"""
File contains some stuff that is used in snapper v3.

author: ulf.schaefer@phe.gov.uk

"""

import logging

from lib.utils import get_closest_threshold
from lib.ClusterStats import ClusterStats
from lib.distances import get_distances

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

def check_merging_needed(cur, distances, new_snad, levels=[0, 5, 10, 25, 50, 100, 250]):
    """
    Checks whether the assignment of the new sample to any of the proposed clusters would require
    the proposed cluster to be merge with any other clusters

    Parameters
    ----------
    cur: obj
        database cursor
    distances: list of tuples
        sorted list of tuples with (sample_id, distance) with closes sample first
        e.g. [(298, 0), (37, 3), (55, 4)]
    new_snad: list
        [t0 or None, t5 or None, t10, t25, t50, t100, t250]

    Returns
    -------
    merges: dict
        {lvl: [list of clusters that need merging - could be >2]}
    """

    merges = {}

    for newcl, lvl in zip(new_snad, levels):
        if newcl == None:
            # we're making a new cluster at this level for the sample so there is naught to merge
            continue
        else:
            # get all samples ids where the sample is <= threshold away
            samples = [s for (s, d) in distances if d <= lvl]
            # get the cluster names for all of these samples
            sql = "SELECT t"+str(lvl)+" FROM sample_clusters WHERE fk_sample_id IN %s"
            cur.execute(sql, (tuple(samples), ))
            # use clustre name as the key in dict to test if they are all the same
            clusters = {r['t'+str(lvl)]: None for r in cur.fetchall()}
            # if there is only one, it means that all the samples <= t from the new samples are in
            # the same cluster => no merge
            # if there is more than one, it measn that the new sample coule go into two different
            # levels clusters at this level and the two clustes need to be merged
            if len(clusters.keys()) > 1:
                merges[lvl] = clusters.keys()

    return merges

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
        {lvl: [list of clusters that need merging - could be >2]}

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

        print "lvl", lvl
        print "clu", clu

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

        print "current_mems1", current_mems

        merge_w_cluster = []
        # if we need to merge
        if merges.has_key(lvl):

            # if there is a merge, we first need to calculate the stats for the newly created merged cluster
            logging.warning("Merge required at level %s between clusters %s. z-score will be checked for the new cluster resulting from this merge!", lvl, str(merges[lvl]))
            oStats, current_mems = get_stats_for_merge(cur, lvl, merges)

        else:
            if nof_mems == 1:
                logging.info("Cluster %s at level %s has only one member. Skipping zscore check.", clu, t_lvl)
                continue

            # if there is no merge we can just take the stats as stored in the db
            oStats = ClusterStats(members=nof_mems, stddev=row['stddev'], mean=row['mean_pwise_dist'])

        print "current_mems2", current_mems

        # get the mean distance of all current members to the new member
        all_dist_to_new_mem = [d for (s, d) in distances if s in current_mems]
        avg_dis = sum(all_dist_to_new_mem) / float(len(current_mems))

        # add a new member to the cluster and update stats
        oStats.add_member(all_dist_to_new_mem)

        # calculate zscore
        mean_all_dist_in_c = oStats.mean_pw_dist
        zscr = (mean_all_dist_in_c - avg_dis) / oStats.stddev_pw_dist

        mess = "z-score of new sample to cluster %s on level %s: %s" % (clu, t_lvl, zscr)
        logging.debug(mess)

        if zscr <= -1.75:
            fail = True
            info.append(mess)

        # do this for all members of the cluster
        nof_mems = len(current_mems)
        for c_mem in current_mems:

            print "c_mem:", c_mem, "t_lvl:", t_lvl

            # get the mean distance of this sample to all other samples in the cluster (w/o the one to be added)
            old_medis = None
            if merges.has_key(lvl) == True:
                # there was a merge, so we can't use what's in the db
                old_medis = get_mean_distance_for_merged_cluster(cur, c_mem, current_mems)
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

            print "old_medis", old_medis
            print "nof_mems", nof_mems
            print "new_dist", new_dist

            # nof_mems is the number of members w/o the sample that we want to add to the cluster
            # update the mean distance with the new distance and calculate the z-score
            new_medis = ((old_medis * (nof_mems - 1)) + new_dist ) / float(nof_mems)

            print "new_medis", new_medis

            zscr = (mean_all_dist_in_c - new_medis) / oStats.stddev_pw_dist

            mess = "z-score of sample %s to cluster %s on level %s incl new member: %s" % (c_mem, clu, t_lvl, zscr)
            # logging.debug(mess)

            if zscr <= -1.75:
                fail = True
                info.append(mess)

    return fail, info

# --------------------------------------------------------------------------------------------------

def get_mean_distance_for_merged_cluster(cur, samid, mems):
    """
    Get the mean distance of a sample (samid) to all samples in the mems list.

    Parameters
    ----------
    cur: obj
        database cursor
    samid: int
        sample_id
    mems: list of int
        all members if this cluster
    Returns
    -------
    m: float
        mean distance
    """

    m = None
    assert samid in mems
    others = [x for x in mems if x != samid]
    dists = get_distances(cur, samid, others)
    d = [d for (s, d) in dists]
    assert len(d) == len(others)
    m = sum(d) / float(len(d))
    return m

# --------------------------------------------------------------------------------------------------

def get_stats_for_merge(cur, lvl, merges):
    """
    Get a stats object for two (or more) clusters after they have been merged:
    either: get the biggest cluster and get the stats from the database
            the add on emember at a time from the other cluster(s)
    or: (if we're merging clusters with only one member) get all pw distances
        in the merged cluster and create stats object with that

    Parameters
    ----------
    cur: obj
        database cursor
    lvl: int
        cluster level (e.g. 10)
    merges: dict
        {lvl: [list of clusters that need merging - could be >2]}

    Returns
    -------
    oStats: obj
        ClusterStats object
    current_mems: list of ints
        list of members of the merged cluster
    None if problem
    """

    # get the members fir each of the clusters to merge and put them in dict
    # members[clu_id] = [list of sample ids]
    clu_to_merge = merges[lvl]
    members = {}
    t_lvl = 't%i' % lvl
    for ctm in clu_to_merge:
        sql = "SELECT fk_sample_id FROM sample_clusters WHERE " + t_lvl + "=%s"
        cur.execute(sql, (ctm, ))
        rows = cur.fetchall()
        ctm_mems = [r['fk_sample_id'] for r in rows]
        members[ctm] = ctm_mems

    # this now has sample_id of the largest cluster first
    clu_to_merge = sorted(members, key=lambda k: len(members[k]), reverse=True)

    # get stats for the biggest cluster from database
    sql = "SELECT nof_members, nof_pairwise_dists, mean_pwise_dist, stddev FROM cluster_stats WHERE cluster_level=%s AND cluster_name=%s"
    cur.execute(sql, (t_lvl, clu_to_merge[0], ))
    if cur.rowcount != 1:
        return None, None
    row = cur.fetchone()

    if row['nof_members'] > 1:

        # make a stats obj from the stuff for the biggest cluster from the db
        oStats = ClusterStats(members=row['nof_members'], stddev=row['stddev'], mean=row['mean_pwise_dist'])

        # get members of biggest cluster
        current_mems = members[clu_to_merge[0]]
        # get all the samples that need adding to it to facilitate the merge
        new_members = []
        for ctm in clu_to_merge[1:]:
            new_members += members[ctm]

        # get all distances for new members and uppdate stats obj iteratively
        for nm in new_members:
            dists = get_distances(cur, nm, current_mems)
            all_dists_to_new_member = [d for (s, d) in dists]
            oStats.add_member(all_dists_to_new_member)
            current_mems.append(nm)

    else:
        # if the biggest cluster has only one member, get all members of all clusters to be merged
        # and get all pw distances - shouldn't be many

        # make a flat list out of the values in members which are lists
        current_mems = [item for sublist in members.values() for item in sublist]

        all_pw_dists = get_all_pw_dists(cur, samids)
        oStats = ClusterStats(members=len(current_mems), dists=all_pw_dists)

    return oStats, current_mems

# --------------------------------------------------------------------------------------------------
