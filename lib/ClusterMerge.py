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
        self.final_name = None
        self.final_members = None
        self.stats = None
        self.member_stats = None

    # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

    def __str__(self):
        '''
        overwrite __str__() to tell it how to print itself
        '''
        return "level: %s, clusters: %s, final name: %s, nof final members: %s" \
                %(self.t_level,
                  self.org_clusters,
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

    def update_tables(self, cur):
        '''
        Update the entries in the tables for this merge.

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

        sql = "UPDATE sample_clusters SET "+self.t_level+"=%s WHERE "+self.t_level+" IN %s"
        cur.execute(sql, (self.final_name, tuple(self.org_clusters), ))

        for fm in self.final_members:
            sql = "UPDATE sample_clusters SET "+self.t_level+"_mean=%s WHERE fk_sample_id=%s"
            cur.execute(sql, (self.member_stats[fm], fm, ))

        return 0

    # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

# --------------------------------------------------------------------------------------------------
