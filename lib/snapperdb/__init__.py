"""
File contains some stuff that is used in snapper v3.

author: ulf.schaefer@phe.gov.uk

"""

import logging

from lib.utils import get_closest_threshold
from lib.ClusterStats import ClusterStats
from lib.distances import get_distances
from lib.merging import get_stats_for_merge, get_mean_distance_for_merged_cluster

# --------------------------------------------------------------------------------------------------

def get_sample_id(cur, name):
    """
    Get the pk_id for a given sample name from the db.
    Indicate error if not unique or not found.

    Parameters
    ----------
    cur: obj
        database cursor
    name: str
        sample name

    Returns
    -------
    sample_id: int
        actual id (>0) if successful
        -1 if fail
    """

    # get sample_id name from database
    sql = "SELECT pk_id FROM samples WHERE sample_name=%s"
    cur.execute(sql, (name, ))
    if cur.rowcount < 1:
        logging.error("A sample with name %s doesn't exist in the database.", name)
        return -1
    elif cur.rowcount > 1:
        logging.error("More than one sample with name %s exist in the database.", name)
        return -1
    else:
        pass
    sample_id = cur.fetchone()[0]

    return sample_id


# --------------------------------------------------------------------------------------------------

def get_closest_samples(cur, distances):
    """
    Get info about the immediate neighbourhood of this sample, i.e. all samples that are closer
    or equally close to the

    Parameters
    ----------
    cur: obj
        database cursor
    distances: list of tuples
        sorted list of tuples with (sample_id, distance) with closes sample first
        e.g. [(298, 0), (37, 3), (55, 4)]

    Returns
    -------
    no name: dict
        {'closest_distance': int,
         'nearest_t': int,
         'closest_sample': int,
         'closest_snad': [list of 7 ints]}
    """

    closest_distance = distances[0][1]
    closest_sample = distances[0][0]

    # returns None if sample is >250 away from closest neighbour
    nearest_t = get_closest_threshold(closest_distance)

    sql = "SELECT t0, t5, t10, t25, t50, t100, t250 FROM sample_clusters WHERE fk_sample_id=%s"
    cur.execute(sql, (closest_sample, ))
    row = cur.fetchone()
    closest_snad = [row['t0'], row['t5'], row['t10'], row['t25'], row['t50'], row['t100'], row['t250']]

    return {'closest_distance': closest_distance,
            'nearest_t': nearest_t,
            'closest_sample': closest_sample,
            'closest_snad': closest_snad}

# --------------------------------------------------------------------------------------------------

def get_new_snp_address(nbhood, levels=[0, 5, 10, 25, 50, 100, 250]):
    """
    Get the proposed new SNP address for a sample based on it's neighbourhood.

    Parameters
    ----------
    nbhood: dist
        {'closest_distance': int,
         'nearest_t': int,
         'closest_sample': int,
         'closest_snad': [list of 7 ints]}

    Returns
    -------
    snad: list
        [t0 or None, t5 or None, t10, t25, t50, t100, t250]
    """

    # get the snpaddr from the closest sample
    snad = nbhood['closest_snad']
    x = 0
    try:
        # overwrite closest sample snpaddr with None until within the threshold
        while nbhood['closest_distance'] > levels[x]:
            snad[x] = None
            x += 1
    except IndexError:
        # this happens when the closest sample is >250 away
        pass
    return snad

# --------------------------------------------------------------------------------------------------

def check_zscores(cur, distances, new_snad, nbhood, merges, levels=[0, 5, 10, 25, 50, 100, 250]):
    """
    Check the zscores of putting a new sample in the clusters proposed, considering merges.

    Parameters
    ----------
    cur: obj
        database cursor
    distances: list of tuples
        sorted list of tuples with (sample_id, distance) with closes sample first
        e.g. [(298, 0), (37, 3), (55, 4)]
    new_snad: list
        [t0 or None, t5 or None, t10, t25, t50, t100, t250]
    nbhood: dict
        {'closest_distance': int,
         'nearest_t': int,
         'closest_sample': int,
         'closest_snad': [list of 7 ints]}
    merges: dict
        {lvl: ClusterMerge object}

    Returns
    -------
    fail: boolean
        False if all successful, True if failure
    info: list
        list of str with info what failed
    None if there is a problem
    """

    fail = False
    info = []

    for (clu, lvl) in zip(new_snad, levels):

        # new cluster at this level, no check required
        if clu == None:
            logging.info("New cluster at t%s level, no zscore check required.", lvl)
            continue

        # get existing stats for this cluster and create ClusterStats obj
        t_lvl = 't%i' % lvl
        sql = "SELECT nof_members, nof_pairwise_dists, mean_pwise_dist, stddev FROM cluster_stats WHERE cluster_level=%s AND cluster_name=%s"
        cur.execute(sql, (t_lvl, clu, ))
        if cur.rowcount != 1:
            return None, None
        row = cur.fetchone()
        nof_mems = row['nof_members']

        # get all current members of this cluster
        sql = "SELECT s.pk_id as samid FROM samples s, sample_clusters c WHERE c."+t_lvl+"=%s AND s.pk_id=c.fk_sample_id AND s.ignore_zscore IS FALSE"
        cur.execute(sql, (clu, ))
        assert cur.rowcount == nof_mems
        rows = cur.fetchall()
        current_mems = [r['samid'] for r in rows]

        oStats = None
        # if we need to merge
        if merges.has_key(lvl):

            # if there is a merge, we first need to calculate the stats for the newly created merged cluster
            logging.warning("Merge required at level %s between clusters %s. z-score will be checked for the new cluster resulting from this merge!", lvl, str(merges[lvl]))
            current_mems = get_stats_for_merge(cur, merges[lvl])
            # Create a new ClusterStats object from the values in the ClusterMerge object.
            # This is because we're adding a member to it below to calculate the zscores.
            # But we don't want that to happend to the stats object in the merge.
            oStats = ClusterStats(members=merges[lvl].stats.members, stddev=merges[lvl].stats.stddev_pw_dist, mean=merges[lvl].stats.mean_pw_dist)

        else:
            if nof_mems == 1:
                logging.info("Cluster %s at level %s has only one member. Skipping zscore check.", clu, t_lvl)
                continue

            # if there is no merge we can just take the stats as stored in the db
            oStats = ClusterStats(members=nof_mems, stddev=row['stddev'], mean=row['mean_pwise_dist'])

        # get the mean distance of all current members to the new member
        all_dist_to_new_mem = [d for (s, d) in distances if s in current_mems]
        avg_dis = sum(all_dist_to_new_mem) / float(len(current_mems))

        # add a new member to the cluster and update stats
        oStats.add_member(all_dist_to_new_mem)

        if oStats.mean_pw_dist <= 0.0:
            logging.info("All distances in cluster %s on level %s are 0. Skipping zscore check.", clu, t_lvl)
            continue

        if oStats.stddev_pw_dist <= 0.0:
            logging.info("All distances in cluster %s on level %s are identical. Skipping zscore check.", clu, t_lvl)
            continue

        # calculate zscore
        mean_all_dist_in_c = oStats.mean_pw_dist
        zscr = (avg_dis - mean_all_dist_in_c) / oStats.stddev_pw_dist

        mess = "z-score of new sample to cluster %s on level %s: %s" % (clu, t_lvl, zscr)
        logging.debug(mess)

        if zscr <= -1.75:
            fail = True
            info.append(mess)

        # do this for all members of the cluster
        nof_mems = len(current_mems)
        merge_per_sample_stats = {}
        for c_mem in current_mems:

            # get the mean distance of this sample to all other samples in the cluster (w/o the one to be added)
            old_medis = None
            if merges.has_key(lvl) == True:
                # there was a merge, so we can't use what's in the db
                old_medis = get_mean_distance_for_merged_cluster(cur, c_mem, current_mems)
                merge_per_sample_stats[c_mem] = old_medis
            else:
                # if there was no merge, get the mean distance of this member to all other members
                # (excluding the one to be added) from the database
                sql = "SELECT "+t_lvl+"_mean FROM sample_clusters WHERE fk_sample_id=%s"
                cur.execute(sql, (c_mem, ))
                if cur.rowcount != 1:
                    return None, None
                row = cur.fetchone()
                old_medis = row[t_lvl + '_mean']

            # get the distance of this member to the sample that we want to add to the cluster
            new_dist = [d for (s, d) in distances if s == c_mem][0]

            # nof_mems is the number of members w/o the sample that we want to add to the cluster
            # update the mean distance with the new distance and calculate the z-score
            new_medis = ((old_medis * (nof_mems - 1)) + new_dist ) / float(nof_mems)

            zscr = (new_medis - mean_all_dist_in_c) / oStats.stddev_pw_dist

            mess = "z-score of sample %s to cluster %s on level %s incl new member: %s" % (c_mem, clu, t_lvl, zscr)
            logging.debug(mess)

            if zscr <= -1.0:
                fail = True
                info.append(mess)

        # if there was a merge we want to remember that we already calculated all this stuff
        if merges.has_key(lvl) == True:
            merges[lvl].member_stats = merge_per_sample_stats

    return fail, info

# --------------------------------------------------------------------------------------------------

def check_duplicate_clustering(cur, sample_id):
    """
    Check whether the sample was already clustered by looking in the
    sample_clusters table. It's bad to sample clusters twice.

    Parameters
    ----------
    cur: obj
        database cursor
    sample_id: int
        pk_id in samples table

    Returns
    -------
    rc: int
        0 if not clustered yet, else <0
    snad: str
        snp address if clustered already, else None

    """

    sql = "SELECT t250, t100, t50, t25, t10, t5, t0 FROM sample_clusters WHERE fk_sample_id=%s"
    cur.execute(sql, (sample_id, ))
    if cur.rowcount == 0:
        rc = 0
        snad = None
    else:
        rc = -1
        snad = '-'.join([str(x) for x in cur.fetchone()])

    return rc, snad

# --------------------------------------------------------------------------------------------------
