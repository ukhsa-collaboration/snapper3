import logging

import argparse
from argparse import RawTextHelpFormatter
import psycopg2
from psycopg2.extras import DictCursor

from lib.distances import get_distances, get_all_pw_dists
from lib.ClusterStats import ClusterStats
from lib.utils import get_all_cluster_members

from datetime import datetime

# --------------------------------------------------------------------------------------------------

def get_desc():
    """
    Get the description of this module
    Parameters
    ----------
    no inputs
    Returns
    -------
    no name: str
        a string containing the description
    """

    return r'''This removes a samples from the database, which is not a trivial thing to do because
cluster stats and stats for other samples need to be updated.

Also: The integrity of all clusters the sampe is in needs to be checked for the potential need to
split them, and if necessary clusters need splitting.

__WARNING__: This may require the calculation of A LOT of distances and MAY TAKE A LONG TIME. Don't
use this unless you have to. You should have thought about it before you put the sample into the
database in the first place.'''

# --------------------------------------------------------------------------------------------------

def get_args():
    """
    Parge arguments
    Parameters
    ----------
    no inputs
    Returns
    -------
    args: obj
        arguments object
    """

    args = argparse.ArgumentParser(description=get_desc(), formatter_class=RawTextHelpFormatter)

    args.add_argument("--connstring",
                      "-c",
                      type=str,
                      metavar="CONNECTION",
                      required=True,
                      dest="db",
                      help="Connection string for db. REQUIRED")

    args.add_argument("--sample",
                      "-s",
                      type=str,
                      required=True,
                      dest="sample",
                      help="Name of sample to remove. REQUIRED.")

    grp = args.add_mutually_exclusive_group()
    grp.add_argument("--just-ignore",
                     action='store_true',
                     help="""Sample and variant information will be retained in database,
but clustering information will be removed. Ignore_sample will
be set to TRUE. [DEFAULT: Remove everything. Sample can be added
and clustered again later.]""")

    grp.add_argument("--known-outlier",
                     action='store_true',
                     help="""Sample, variant, and clustering information will be retained in database,
but cluster stats will be reverted. ignore_zscore will
be set to TRUE. [DEFAULT: Remove everything. Sample can be added
and clustered again later.]""")


    return args

# --------------------------------------------------------------------------------------------------

def main(args):
    '''
    Main funtion
    Parameters
    ----------
    no inputs
    Returns
    -------
    0
    Creates all logs and result files
    '''

    logging.debug("Args received: %s", str(args))

    try:
        # open db
        conn = psycopg2.connect(args['db'])
        cur = conn.cursor(cursor_factory=psycopg2.extras.DictCursor)

        # get the pkid of this sample
        sql = "SELECT pk_id, ignore_sample, ignore_zscore FROM samples WHERE sample_name=%s"
        cur.execute(sql, (args['sample'], ))
        if cur.rowcount < 1:
            logging.error("Sample '%s' not found in database.", args['sample'])
            return 1
        row = cur.fetchone()

        sample_id = row['pk_id']
        # get status of ignore flags
        igno_flag = row['ignore_sample']
        zscr_flag = row['ignore_zscore']

        # this sample is already ignored
        if igno_flag == True:
            logging.info("This sample is already ignored.")
            if args['just_ignore'] == True:
                logging.info("You chose not to remove it completely.")
                return 0
            elif args['known_outlier'] == True:
                logging.error("This sample is not clustered, so can't be made an outlier.")
                return 1
            else:
                drop_sample(cur, sample_id)
                conn.commit()
                return 0

        # get the snp address
        sql = "SELECT t0, t5, t10, t25, t50, t100, t250 FROM sample_clusters WHERE fk_sample_id=%s"
        cur.execute(sql, (sample_id, ))
        if cur.rowcount > 1:
            logging.error("There is not exactly one entry in sample clusters for this samples. :-(")
            return 1
        elif cur.rowcount == 0:
            logging.info("Sample has never been clustered.")
            if args['just_ignore'] == True:
                sql = "UPDATE samples SET ignore_sample=True WHERE pk_id=%s"
                cur.execute(sql, (sample_id, ))
                logging.info("Sample is now ignored.")
            elif args['known_outlier'] == True:
                logging.error("A sample not clustered cannot be a 'known outlier'. Please select to ignore or remove the sample.")
                return 1
            else:
                drop_sample(cur, sample_id)
            conn.commit()
            return 0
        row = cur.fetchone()
        snad = [row['t0'], row['t5'], row['t10'], row['t25'], row['t50'], row['t100'], row['t250']]

        # 'global' container to keep distances to prevent repeated calculation of same distance
        # always keep distances[a][b] = d and distances[b][a] = d both
        distances = {}

        # this sample is a know outlier
        if zscr_flag == True:
            logging.info("This sample is a known outlier.")
            if args['just_ignore'] == True:

                if update_clustering(cur, sample_id, snad, distances, zscr_flag) != 0:
                    logging.error("Error in update clustering.")
                    return 1

                logging.info("Removing the sample from the sample_clusters and sample_history table.")
                sql = "DELETE FROM sample_clusters WHERE fk_sample_id=%s"
                cur.execute(sql, (sample_id, ))
                sql = "DELETE FROM sample_history WHERE fk_sample_id=%s"
                cur.execute(sql, (sample_id, ))
                sql = "UPDATE samples SET (ignore_sample, ignore_zscore) = (True, False) WHERE pk_id=%s"
                cur.execute(sql, (sample_id, ))
                logging.info("Sample is now ignored.")
                conn.commit()
            elif args['known_outlier'] == True:
                logging.info("So there is nothing to do.")
            else:

                if update_clustering(cur, sample_id, snad, distances, zscr_flag) != 0:
                    logging.error("Error in update clustering.")
                    return 1

                logging.info("Removing the sample from the sample_clusters and sample_history table.")
                sql = "DELETE FROM sample_clusters WHERE fk_sample_id=%s"
                cur.execute(sql, (sample_id, ))
                sql = "DELETE FROM sample_history WHERE fk_sample_id=%s"
                cur.execute(sql, (sample_id, ))
                drop_sample(cur, sample_id)
                conn.commit()
            return 0

        # ---------------------------------------------------------
        # if we get to here we know the sample is 'fully clustered'
        # ---------------------------------------------------------

        if args['just_ignore'] == True:

            if update_clustering(cur, sample_id, snad, distances, zscr_flag) != 0:
                logging.error("Error in update clustering.")
                return 1

            # now remove the sample
            logging.info("Removing the sample from the sample_clusters and sample_history table.")
            sql = "DELETE FROM sample_clusters WHERE fk_sample_id=%s"
            cur.execute(sql, (sample_id, ))
            sql = "DELETE FROM sample_history WHERE fk_sample_id=%s"
            cur.execute(sql, (sample_id, ))
            logging.info("You chose not to remove the sample completely.")
            sql = "UPDATE samples SET ignore_sample=True WHERE pk_id=%s"
            cur.execute(sql, (sample_id, ))

        elif args['known_outlier'] == True:
            if make_known_outlier(cur, sample_id, snad, distances) != 0:
                logging.error("Error in update clustering.")
                return 1
            logging.info("Sample successfully turned into a known outlier.")
        else:
            if update_clustering(cur, sample_id, snad, distances, zscr_flag) != 0:
                logging.error("Error in update clustering.")
                return 1

            # now remove the sample
            logging.info("Removing the sample from the sample_clusters and sample_history table.")
            sql = "DELETE FROM sample_clusters WHERE fk_sample_id=%s"
            cur.execute(sql, (sample_id, ))
            sql = "DELETE FROM sample_history WHERE fk_sample_id=%s"
            cur.execute(sql, (sample_id, ))
            drop_sample(cur, sample_id)

        conn.commit()

    except psycopg2.Error as e:
         logging.error("Database reported error: %s" % (str(e)))
         return None
    finally:
        # close all dbs
        cur.close()
        conn.close()

    return 0

# end of main --------------------------------------------------------------------------------------

def make_known_outlier(cur, sample_id, snad, distances):
    """
    Turns a sample from a fully cluster one to a 'known outlier'. This involves now splitting since
    the sample is kept in the cluster. Stats for the clusters this sample is in need to be updated.

    Parameters
    ----------
    cur: obj
        database cursor
    sample_id: int
        id of sample to remove
    snad: list
        [t0, t5, ..., t250]
    distances: dist
        distances[a][b] = d
        distances[b][a] = d

    Returns
    -------
    1 on error, 0 on success

    """

    # get all 250 members and calculate distances to them
    sql = "SELECT c.fk_sample_id AS samid FROM sample_clusters c, samples s WHERE c.t250=%s AND c.fk_sample_id=s.pk_id AND s.ignore_zscore IS FALSE"
    cur.execute(sql, (snad[6], ))
    t250_members = [r['samid'] for r in cur.fetchall()]
    try:
        t250_members.remove(sample_id)
    except ValueError:
        logging.error("Bizzare data inconsistency for sample id %s and t250 cluster %s.", sample_id, snad[6])
        return  1

    logging.info("Calculating %i distances to update stats.", len(t250_members))
    # get the distances and remember them in the dict
    # we need them anyway, so we store them in the dict
    _ = get_distances_from_memory(cur, distances, sample_id, t250_members)

    for clu, lvl in zip(snad, [0, 5, 10, 25, 50, 100, 250]):
        t_lvl = "t%s" % (lvl)

        # get stats for this cluster
        sql = "SELECT nof_members, nof_pairwise_dists, mean_pwise_dist, stddev FROM cluster_stats WHERE cluster_level=%s AND cluster_name=%s"
        cur.execute(sql, (t_lvl, clu, ))
        if cur.rowcount == 0:
            logging.error("Cluster stats for level %s and cluster %s not found.", t_lvl, clu)
            return 1
        row = cur.fetchone()

        # when deleting the last member of this cluster there is no need to update anything
        if row['nof_members'] <= 1:
            logging.debug("This is the last member of cluster %s on level %s. Cluster now has 0 members.", clu, t_lvl)
            sql = "UPDATE cluster_stats SET (nof_members, nof_pairwise_dists, mean_pwise_dist, stddev) = (0, 0, null, null) WHERE cluster_level=%s AND cluster_name=%s"
            cur.execute(sql, (t_lvl, clu, ))
            continue

        # create cluster stats object from the information in the database
        oStats = ClusterStats(members=row['nof_members'], stddev=row['stddev'], mean=row['mean_pwise_dist'])

        # get other members of this cluster, we know it must be at least one
        sql = "SELECT c.fk_sample_id AS samid FROM sample_clusters c, samples s WHERE c."+t_lvl+"=%s AND c.fk_sample_id=s.pk_id AND s.ignore_zscore IS FALSE"
        cur.execute(sql, (clu, ))
        members = [r['samid'] for r in cur.fetchall()]
        try:
            members.remove(sample_id)
        except ValueError:
            logging.error("Bizzare data inconsistency for sample id %s and %s cluster %s.", sample_id, t_lvl, clu)
            return 1

        this_di = [distances[sample_id][m] for m in members]
        assert oStats.members == (len(this_di) + 1)
        oStats.remove_member(this_di)

        # update the cluster stats in the database with the info from the object
        sql = "UPDATE cluster_stats SET (nof_members, nof_pairwise_dists, mean_pwise_dist, stddev) = (%s, %s, %s, %s) WHERE cluster_level=%s AND cluster_name=%s"
        cur.execute(sql, (oStats.members, oStats.nof_pw_dists, oStats.mean_pw_dist, oStats.stddev_pw_dist, t_lvl, clu, ))

    # for known outliers the mean distances to all other samples in the clusters are not kept
    sql = "UPDATE sample_clusters SET (t0_mean, t5_mean, t10_mean, t25_mean, t50_mean, t100_mean, t250_mean) = (null, null, null, null, null, null, null) WHERE fk_sample_id=%s"
    cur.execute(sql, (sample_id, ))

    # set the flag
    sql = "UPDATE samples SET ignore_zscore=True WHERE pk_id=%s"
    cur.execute(sql, (sample_id, ))

    return 0

# --------------------------------------------------------------------------------------------------

def update_clustering(cur, sample_id, snad, distances, zscr_flag):
    """
    Update the sample clustering. Check for splits.

    Parameters
    ----------
    cur: obj
        database cursor
    sample_id: int
        id of sample to remove
    snad: list
        [t0, t5, ..., t250]
    distances: dist
        distances[a][b] = d
        distances[b][a] = d
    zscr_flag: boolean
        True: sample is 'known outlier', False: sample is 'fully clustered'

    Returns
    -------
    1 on error, 0 on success

    """

    logging.info("Checking if any clusters need splitting after sample removal.")
    splits = check_cluster_integrity(cur, sample_id, snad, distances)
    logging.info("Clusters that need splitting after sample removal: %s", splits)

    # splits structure:
    # splits[level] = [(c, a, b), ...] <- a and b are samples that are no longer connected the
    # cluster c after the removee has been removed
    # splits[level] = None <- at this level the cluster is fine

    # if the sample was previously included in the stats or there are any splits on any level
    # we need to (do the splitting and) update the stats
    if zscr_flag == False or any([x != None for x in splits.values()]):
        # get all 250 members and calculate distances to them
        sql = "SELECT c.fk_sample_id AS samid FROM sample_clusters c, samples s WHERE c.t250=%s AND c.fk_sample_id=s.pk_id AND s.ignore_zscore IS FALSE"
        cur.execute(sql, (snad[6], ))
        t250_members = [r['samid'] for r in cur.fetchall()]
        try:
            logging.debug("Removing %s from list: %s", sample_id, t250_members)
            t250_members.remove(sample_id)
        except ValueError:
            if zscr_flag == True:
                logging.debug("Could not find %s as a member of t250 cluster %s, but it's OK because it's a known outlier.", sample_id, snad[6])
            else:
                logging.error("Bizzare data inconsistency for sample id %s and t250 cluster %s.", sample_id, snad[6])
                return 1

        logging.info("Calculating %i distances to update stats.", len(t250_members))
        # get the distances and remember them in the dict
        # we need them anyway, so we store them in the dict
        _ = get_distances_from_memory(cur, distances, sample_id, t250_members)

        # update all stats in all clusters on all levels
        for clu, lvl in zip(snad, [0, 5, 10, 25, 50, 100, 250]):
            if update_cluster_stats_post_removal(cur, sample_id, clu, lvl, distances, splits[lvl], zscr_flag) == None:
                logging.error("Problem with updating cluster stats.")
                return 1
    else:
        logging.info("No stats update required for this sample.")

    return 0

# --------------------------------------------------------------------------------------------------

def split_clusters(cur, sample_id, problems, lvl, distances):
    """
    Split cluster.

    Parameters
    ----------
    cur: obj
        database cursor
    sample_id: int
        id of sample to remove
    problems: list of tuples
        [(c, a, b), ...] <- no longer connected pairs in cluster c
    distances: dist
        distances[a][b] = d
        distances[b][a] = d

    Returns
    -------
    groups: dict
        groups[a] = [list of others]
        groups[b] = [list of others]

    """

    groups = {}
    visited = set()
    for (c, a, b) in problems:
        mems = get_all_cluster_members(cur, c, 't'+str(lvl))
        for node in [a,b]:
            # get a set of all samples that are already in one of the groups
            _ = [visited.update(x) for x in groups.values()]
            # if we already expanded from that node or went past it in a previous group, don't go
            if groups.has_key(node) or (node in visited):
                continue
            else:
                groups[node] = expand_from_node(cur, node, c, lvl, distances, sample_id)

        # if the combined langth of all groups covers the whole cluster (without the removee), we're done
        if sum([len(x) for x in groups.values()]) == len(mems) - 1:
            break
        # else there must be another 'broken link'

    return groups

# --------------------------------------------------------------------------------------------------

def update_cluster_stats_post_removal(cur, sid, clu, lvl, distances, split, zscr_flag):
    """
    Update the cluster stats and the sample stats for removing the sample from the cluster.

    Parameters
    ----------
    cur: obj
        database cursor
    sid: int
        pk id of sample to remove
    clu: int
        cluster id to remove from
    lvl: int
        cluster level of removal
    distances: list of tuples
        [(samid, distance), (samid, distance), (samid, distance), ...]
    split: list of tuples
        [(c, a, b), ...] <- no longer connected pair in cluster c
    Returns
    -------
    0 if fine
    None if fail
    """

    t_lvl = "t%i" % (lvl)

    logging.info("Updating stats for cluster %s on level %s.", clu, t_lvl)

    # get the cluster stats from the database
    sql = "SELECT nof_members, mean_pwise_dist, stddev FROM cluster_stats WHERE cluster_level=%s AND cluster_name=%s"
    cur.execute(sql, (t_lvl, clu, ))
    if cur.rowcount == 0:
         if zscr_flag == True:
            logging.info("Sample is a known outlier and the only member of level %s cluster %s. So there are no stats to update.", t_lvl, clu)
            return 0
         else:
            logging.error("Cluster stats for level %s and cluster %s not found.", t_lvl, clu)
            return None
    row = cur.fetchone()

    # when deleting the last member of this cluster there is no need to update anything, just get rid of it
    if row['nof_members'] <= 1:
        logging.debug("This is the last member of cluster %s on level %s. Deleting cluster stats.", clu, t_lvl)
        sql = "DELETE FROM cluster_stats WHERE cluster_level=%s AND cluster_name=%s"
        cur.execute(sql, (t_lvl, clu, ))
        return 0

    # create cluster stats object from the information in the database
    oStats = ClusterStats(members=row['nof_members'], stddev=row['stddev'], mean=row['mean_pwise_dist'])

    # get other members of this cluster, we know it must be at least one
    sql = "SELECT c.fk_sample_id AS samid FROM sample_clusters c, samples s WHERE c."+t_lvl+"=%s AND c.fk_sample_id=s.pk_id AND s.ignore_zscore IS FALSE"
    cur.execute(sql, (clu, ))
    members = [r['samid'] for r in cur.fetchall()]
    logging.debug("Got the follwoing members: %s", members)
    try:
        logging.debug("Removing %s from list: %s", sid, members)
        members.remove(sid)
    except ValueError:
        if zscr_flag == True:
            logging.debug("Could not find %s as a member of %s cluster %s, but it's OK because it's a known outlier.",
                          sid, t_lvl, clu)
        else:
            logging.error("Bizzare data inconsistency for sample id %s and %s cluster %s.", sid, t_lvl, clu)
            return None

    # if the was previously ignore do not remove it from stats object, because it was never considered when calculating the stats
    if zscr_flag == False:
        # get all distances from the sample to be removed to all other members of the cluster
        # and update the stats object with this information
        this_di = [distances[sid][m] for m in members]
        assert oStats.members == (len(this_di) + 1)
        oStats.remove_member(this_di)
        # remember which members have been removed from the stats object
        removed_members = [sid]
    else:
        removed_members = []

    # if tere is a split on this level
    if split != None:
        logging.info("Cluster %s need to be split.", clu)
        groups = split_clusters(cur, sid, split, lvl, distances)
        logging.debug("It will be split into these subclusters: %s", groups)

        # groups[a] = [1,2,3]
        # groups[b] = [4,5,6]

        knwntlrs = set()

        # put the largest subcluster at the front of the list of subclusters
        group_lists = sorted(groups.values(), key=len, reverse=True)
        logging.debug("These are the group lists: %s", group_lists)
        # for the largest group
        for grli in group_lists[1:]:
            # for all members of this group
            logging.debug("Current group list: %s", grli)
            for m in grli:
                # remove from members, from stats object and remember that you removed it in that list
                logging.debug("Removing %s from %s", m, members)
                if m in members:
                    members.remove(m)
                    this_di = [d for (s, d) in get_distances_from_memory(cur, distances, m, members)]
                    assert oStats.members == (len(this_di) + 1)
                    # i.e. turn the oStats object into the stats object for the largest cluster after the split
                    oStats.remove_member(this_di)
                    removed_members.append(m)
                else:
                    knwntlrs.add(m)
                    logging.debug("Could not remove %s from members, but it's probably a kown outlier.")

        # for the other subclustrs
        for grli in group_lists[1:]:
            # make a new stats object based on the list of members and all pw distances between them
            # remove known outliers previously encountered from consideration
            grli = list(set(grli).difference(knwntlrs))
            all_pw_grdi = get_all_pw_dists(cur, grli)
            oStatsTwo = ClusterStats(members=len(grli), dists=all_pw_grdi)
            sql = "SELECT max("+t_lvl+") AS m FROM sample_clusters"
            cur.execute(sql)
            row = cur.fetchone()
            new_clu_name = row['m'] + 1
            sql = "INSERT INTO cluster_stats (cluster_level, cluster_name, nof_members, nof_pairwise_dists, mean_pwise_dist, stddev) VALUES (%s, %s, %s, %s, %s, %s)"
            cur.execute(sql, (t_lvl, new_clu_name, oStatsTwo.members, oStatsTwo.nof_pw_dists, oStatsTwo.mean_pw_dist, oStatsTwo.stddev_pw_dist, ))

            # document the upcoming change in the sample history
            update_sample_history(cur, t_lvl, new_clu_name, grli)

            # put all members of this subcluster in the new cluster in the database
            sql = "UPDATE sample_clusters SET "+t_lvl+"=%s WHERE fk_sample_id IN %s"
            cur.execute(sql, (new_clu_name, tuple(grli), ))

            # calculate the mean distance to all other members from scratch for all members of this newly
            # created subcluster and update in the database
            for nm in grli:
                targets = [x for x in grli if x != nm]
                alldis = [di for (sa, di) in get_distances_from_memory(cur, distances, nm, targets)]
                try:
                    mean = sum(alldis) / float(len(alldis))
                except ZeroDivisionError:
                    mean = None
                sql = "UPDATE sample_clusters SET "+t_lvl+"_mean=%s WHERE fk_sample_id=%s"
                cur.execute(sql, (mean, nm, ))

    # then update the cluster stats in the database with the info from the object
    sql = "UPDATE cluster_stats SET (nof_members, nof_pairwise_dists, mean_pwise_dist, stddev) = (%s, %s, %s, %s) WHERE cluster_level=%s AND cluster_name=%s"
    cur.execute(sql, (oStats.members, oStats.nof_pw_dists, oStats.mean_pw_dist, oStats.stddev_pw_dist, t_lvl, clu, ))

    # for all other members of this cluster
    for mem in members:

        # get the mean distance to all other members
        sql = "SELECT "+t_lvl+"_mean FROM sample_clusters WHERE fk_sample_id=%s"
        cur.execute(sql, (mem, ))
        if cur.rowcount == 0:
            logging.error("Cluster %s not found in sample_clusters table.", mem)
            return None
        row = cur.fetchone()
        p_mean = row[t_lvl+'_mean']
        n_mean = p_mean

        # update this mean by removing one distance and update the database table
        # if there was not split removed_members will have only one member in it
        # else there are more
        # If there was no split on this level and the samples was previously ignored,
        # it will be empty
        for remomem in removed_members:
            x = get_distances_from_memory(cur, distances, mem, [remomem])[0][1]
            try:
                n_mean = ((p_mean * len(members)) - x) / float(len(members) - 1)
            except ZeroDivisionError:
                n_mean = None
            p_mean = n_mean

        sql = "UPDATE sample_clusters SET "+t_lvl+"_mean=%s WHERE fk_sample_id=%s"
        cur.execute(sql, (n_mean, mem, ))

    return 0

# --------------------------------------------------------------------------------------------------

def drop_sample(cur, sid):
    """
    Remove the sample with the id from the variants and samples tables.

    Parameters
    ----------
    cur: obj
        database cursor
    sid: int
        pk id of sample to remove

    Returns
    -------
    0
    """

    logging.info("Removing the sample from the samples and variants table.")

    sql = "DELETE FROM variants WHERE fk_sample_id=%s"
    cur.execute(sql, (sid, ))

    sql = "DELETE FROM samples WHERE pk_id=%s"
    cur.execute(sql, (sid, ))

    return 0

# --------------------------------------------------------------------------------------------------

def check_cluster_integrity(cur, sample_id, snad, distances, levels=[0, 5, 10, 25, 50, 100, 250]):
    """
    Check whether the removal of sample_id from any of its cluster necessitates
    the split of the cluster.

    Parameters
    ----------
    cur: obj
        database cursor
    sample_id: int
        id of sample to remove
    snad: list of 7 int
        snip address
    distances: dist
        distances[a][b] = d
        distances[b][a] = d
    levels: list of 7 int
        better not change this
        [0, 5, 10, 25, 50, 100, 250]

    Returns
    -------
    None if no splits required, else:
    splits: dict
        splits[level] = [(c, a, b), ...] <- no longer connected pair in cluster c
    """

    splits = {}

    for clu, lvl in zip(snad, levels):

        t_lvl = 't%i' % (lvl)

        logging.info("Checking cluster integrity for cluster %s on level %s.", clu, t_lvl)

        # get all other members of the cluster apart from the removee
        mems = get_all_cluster_members(cur, clu, t_lvl)
        mems.remove(sample_id)

        # get distances of the removee to them
        d = get_distances(cur, sample_id, mems)
        connected_mems = []
        for (sa, di) in d:
            # get all samples that are connected to the removee with d <= t
            if di <= lvl:
                connected_mems.append(sa)
            remember_distance(distances, sample_id, sa, di)

        logging.debug("Samples connected via removee: %s", sorted(connected_mems))

        # investigate all pw distances between connected members
        potentially_broken_pairs = []
        for i, a in enumerate(connected_mems):
            for j, b in enumerate(connected_mems):
                if i < j:
                    pwd = None
                    try:
                        pwd = distances[a][b]
                    except KeyError:
                        pwd = get_all_pw_dists(cur, [a, b])[0]
                        remember_distance(distances, a, b, pwd)
                    # if pw distance between the two sampes is bigger than the threshold,
                    # the link between the samples might be broken, unless there is another samples
                    # (or chain of samples) connecting them
                    if pwd > lvl:
                        potentially_broken_pairs.append((a, b))

        # all pairs that were connected through the removee are also directly connected, happy days
        if len(potentially_broken_pairs) == 0:
            splits[lvl] = None
            continue

        logging.debug("Samples potentially no longer connected via removee: %s",
                      potentially_broken_pairs)

        # check if there is another path to get from a to b with only steps <= t
        for a, b in potentially_broken_pairs:
            broken = False
            logging.debug("Checking if there is another way to connect %s and %s.", a, b)
            # list of samples connectable to a (directly or over multiple nodes)
            rel_conn_sams_to_a = [a]
            idx = 0
            # when b in connected the a w're done
            while b not in rel_conn_sams_to_a:
                # pivot is the one currently investigated
                pivot = rel_conn_sams_to_a[idx]
                # get all the members of the current cluster except the pivot
                all_mems_but_pivot = [x for x in mems if x != pivot]
                # get all the distances from the pivot to thpse members
                d = get_distances_from_memory(cur,
                                              distances,
                                              pivot,
                                              all_mems_but_pivot)
                # all new samples that are connectable to the pivot
                # two conditions: a) the sample is connected to the pivot with d<=t
                #         b) we don't have this sample yet in the ones we already know are connected to a
                rel_conn_to_pivot = [sa for (sa, di) in d
                                     if (di <= lvl) and
                                     (sa not in rel_conn_sams_to_a)]
                # there are no new samples connected to the pivot and the last sample has been considered
                # but b is not yet foud to be connected => cluster is broken
                if len(rel_conn_to_pivot) == 0 and pivot == rel_conn_sams_to_a[-1]:
                    broken = True
                    break
                else:
                    # otehr wise add any potential new ones to the list and check the next one
                    rel_conn_sams_to_a += rel_conn_to_pivot
                    idx += 1
            # we need to remember what was broken for updating later
            if broken == True:
                try:
                    splits[lvl].append((clu, a, b))
                except KeyError:
                    splits[lvl] = [(clu, a, b)]
                # go to next broken pair, there might be more than one and
                # we want to know for updating later

        # we checked all pairs and always found b somehow, cluster is fine
        if splits.has_key(lvl) == False:
            splits[lvl] = None

    return splits

# --------------------------------------------------------------------------------------------------

def get_distances_from_memory(cur, distances, a, targets):
    """
    Get all the distances from 'a' to the target list. Check if they are in
    distances before calculating them. Put the newly calculatd into distances.

    Parameters
    ----------
    cur: obj
        database cursor
    distances: dist
        distances[a][b] = d
        distances[b][a] = d
    a: int
        sample id
    targets: list of int
        list of sample ids

    Returns
    -------
    result: list of tuples
        sorted list of tuples with (sample_id, distance) with closes sample first
        e.g. [(298, 0), (37, 3), (55, 4)]

    """

    result = []
    others = []
    for t in targets:
        try:
            result.append((t, distances[a][t]))
        except KeyError:
            others.append(t)

    if len(others) > 0:
        d = get_distances(cur, a, others)
        for (sa, di) in d:
            remember_distance(distances, a, sa, di)
        result += d

    result = sorted(result, key=lambda x: x[1])
    return result

# --------------------------------------------------------------------------------------------------

def remember_distance(distances, a, b, d):
    """
    Store the distance d from a to b in distances dictionary

    Parameters
    ----------
    distances: dist
        distances[a][b] = d
        distances[b][a] = d
    a: int
        sample id
    b: int
        sample id
    d: int
        distance from a to b

    Returns
    -------
    None
    """

    try:
        distances[a][b] = d
    except KeyError:
        distances[a] = {b: d}
    try:
        distances[b][a] = d
    except KeyError:
        distances[b] = {a: d}
    return None

# --------------------------------------------------------------------------------------------------

def expand_from_node(cur, a, c, lvl, distances, sample_id=None):
    """
    From a seed samples a, get all samples that can be connect with <=lvl over multiple nodes.

    Parameters
    ----------
    cur: obj
        database cursor
    a: int
        seed samples id
    c: int
        cluster name
    lvl: int
        cluster level
    distances: dict
        distances[a][b] = d
        distances[b][a] = d
    sample_id: int
        (optional)  sample to remove from cluster

    Returns
    -------
    with_a: list
        all samples connected to a (over multiple nodes)
    """

    logging.info("Expanding from sample %s.", a)
    t_lvl = 't%i' % (lvl)
    with_a = [a]
    mems = get_all_cluster_members(cur, c, t_lvl)
    # optionally we can remove one
    if sample_id != None:
        mems.remove(sample_id)

    idx = 0
    while True:
        pivot = with_a[idx]
        dis = get_distances_from_memory(cur, distances, pivot, mems)
        # get all the new nodes with are connected to the curren pivot and which we have not considered yet
        new_nodes = [sa for (sa, di) in dis if di <= lvl and sa not in with_a]
        # if there are no new ones and we have reached the end of the list, we're done
        if len(new_nodes) == 0 and pivot == with_a[-1]:
            break
        else:
            # else add more nodes to consider
            with_a += new_nodes
        # consider the next one
        idx += 1

    logging.info("Samples connected to sample %s: %s.", a, with_a)

    return with_a

# --------------------------------------------------------------------------------------------------

def update_sample_history(cur, lvl_in, new_clu_name, grli):
    """
    When changing a snp address because the cluster was split, this need to go into sample_history.

    Parameters
    ----------
    cur: obj
        database cursor
    t_lvl:
        level where snp address changed, e.g. 't100'
    new_clu_name:
        the new cluster name at this level
    grli:
        list of sampleids for which this change is effective

    Returns
    -------
    0
    """

    levels = ['t0', 't5', 't10', 't25', 't50', 't100', 't250']

    for samid in grli:

        sql = "SELECT t0, t5, t10, t25, t50, t100, t250 FROM sample_clusters WHERE fk_sample_id=%s"
        cur.execute(sql, (samid, ))
        oldsnad = cur.fetchone()

        # one line conditional dictionary comprehension
        newsnad = {lvl: oldsnad[lvl] if lvl != lvl_in else new_clu_name for lvl in levels}

        sql = "INSERT INTO sample_history (fk_sample_id, t250_old, t100_old, t50_old, t25_old, t10_old, t5_old, t0_old, t250_new, t100_new, t50_new, t25_new, t10_new, t5_new, t0_new, renamed_at) VALUES (%s, %s, %s, %s, %s, %s, %s, %s, %s, %s, %s, %s, %s, %s, %s, %s)"
        cur.execute(sql, (samid,
                          oldsnad['t250'], oldsnad['t100'], oldsnad['t50'], oldsnad['t25'],
                          oldsnad['t10'], oldsnad['t5'], oldsnad['t0'],
                          newsnad['t250'], newsnad['t100'], newsnad['t50'], newsnad['t25'],
                          newsnad['t10'], newsnad['t5'], newsnad['t0'],
                          datetime.now()))

    return 0

# --------------------------------------------------------------------------------------------------

if __name__ == "__main__":
    exit(main(vars(get_args().parse_args())))
