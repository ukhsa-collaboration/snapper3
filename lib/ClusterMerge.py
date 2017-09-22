#!/usr/bin/env python
"""

"""
import sys
import logging
from datetime import datetime

from lib.distances import get_distances

__version__= '0.1'
__date__= '14Jul2017'
__author__ = 'ulf.schaefer@phe.gov.uk'

# --------------------------------------------------------------------------------------------------

class ClusterMerge(object):
    '''
    Class to help with cluster merging.

    '''
    def __init__(self, **kwargs):
        """
        Constructor of the class. Requires two kwarg:
            - level (int)
            - clusters (list of int)
        """

        self.level = int(kwargs['level'])
        self.t_level = 't%i' % self.level
        self.org_clusters = kwargs['clusters']
        self.org_sizes = kwargs['sizes']
        self.final_name = None
        self.final_members = None
        self.stats = None
        self.member_stats = None

    # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

    def __str__(self):
        '''
        overwrite __str__() to tell it how to print itself
        '''
        return "level: %s, clusters: %s, cluster_sizes: %s, final name: %s, nof final members: %s" \
                %(self.t_level,
                  self.org_clusters,
                  self.org_sizes,
                  "tbd" if self.final_name == None else self.final_name,
                  "tbd" if self.final_members == None else len(self.final_members))

    # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

    def calculate_per_member_stats(self, cur):
        '''
        Calculate the mean distance of all members of the newly merged cluster to all other members

        Parameters
        ----------
        self: obj
            that's me
        cur: obj
            database cursor

        Returns
        -------
        0
        '''

        self.member_stats = {}

        logging.info("Calculating mean distance of all members of merging cluster %s on level %s.", self.final_name, self.t_level)

        for fm in self.final_members:
            others = [x for x in self.final_members if x != fm]
            dists = get_distances(cur, fm, others)
            x = [d for (s, d) in dists]
            self.member_stats[fm] = sum(x) / float(len(x))

        return 0

    # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

    def update_tables(self, cur, levels=['t250', 't100', 't50', 't25', 't10', 't5', 't0']):
        '''
        Update the entries in the tables for this merge.

        Parameters
        ----------
        self: obj
            that's me
        cur: obj
            database cursor
        levels: list
            default: ['t250', 't100', 't50', 't25', 't10', 't5', 't0']

        Returns
        -------
        0
        '''

        sql = "UPDATE cluster_stats SET (nof_members, nof_pairwise_dists, mean_pwise_dist, stddev) = (%s, %s, %s, %s) WHERE cluster_level=%s AND cluster_name=%s"
        cur.execute(sql, (self.stats.members,
                          self.stats.nof_pw_dists,
                          self.stats.mean_pw_dist,
                          self.stats.stddev_pw_dist,
                          self.t_level,
                          self.final_name, ))

        clu_to_del = set(self.org_clusters).difference(set([self.final_name]))

        sql = "DELETE FROM cluster_stats WHERE cluster_level=%s AND cluster_name IN %s"
        cur.execute(sql, (self.t_level, tuple(clu_to_del), ))

        logging.warning("The clusters %s on level %s have been MERGED into cluster %s and have been DELETED.", str(list(clu_to_del)), self.t_level, self.final_name)

        nownow = datetime.now()
        for source in clu_to_del:
            sql = "INSERT INTO merge_log (cluster_level, source_cluster, target_cluster, time_of_merge) VALUES (%s, %s, %s, %s)"
            cur.execute(sql, (self.t_level, source, self.final_name, nownow, ))

        # write to the log which samples get changed from what to what
        sql = "SELECT s.pk_id, s.sample_name, c.t0, c.t5, c.t10, c.t25, c.t50, c.t100, c.t250 FROM samples s, sample_clusters c WHERE s.pk_id=c.fk_sample_id AND c."+self.t_level+" IN %s"
        cur.execute(sql , (tuple([x for x in self.org_clusters if x != self.final_name]), ))
        rows = cur.fetchall()
        for r in rows:
            logging.warning("Clustering for sample %s will be changed from %s to %s",
                            r['sample_name'],
                            '-'.join([str(r[x]) for x in levels]),
                            '-'.join([str(self.final_name) if x == self.t_level else str(r[x]) for x in levels]))

            # also write this to a table
            sql = "INSERT INTO sample_history (fk_sample_id, t250_old, t100_old, t50_old, t25_old, t10_old, t5_old, t0_old, t250_new, t100_new, t50_new, t25_new, t10_new, t5_new, t0_new, renamed_at) VALUES (%s, %s, %s, %s, %s, %s, %s, %s, %s, %s, %s, %s, %s, %s, %s, %s)"
            cur.execute(sql ,(r['pk_id'], \
                              r['t250'], r['t100'], r['t50'], r['t25'], r['t10'], r['t5'], r['t0'], \
                              self.final_name if self.t_level == 't250' else r['t250'], \
                              self.final_name if self.t_level == 't100' else r['t100'], \
                              self.final_name if self.t_level == 't50' else r['t50'], \
                              self.final_name if self.t_level == 't25' else r['t25'], \
                              self.final_name if self.t_level == 't10' else r['t10'], \
                              self.final_name if self.t_level == 't5' else r['t5'], \
                              self.final_name if self.t_level == 't0' else r['t0'], \
                              nownow, ))

        sql = "UPDATE sample_clusters SET "+self.t_level+"=%s WHERE "+self.t_level+" IN %s"
        cur.execute(sql, (self.final_name, tuple(self.org_clusters), ))

        for fm in self.final_members:
            sql = "UPDATE sample_clusters SET "+self.t_level+"_mean=%s WHERE fk_sample_id=%s"
            cur.execute(sql, (self.member_stats[fm], fm, ))

        return 0

    # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

# --------------------------------------------------------------------------------------------------
