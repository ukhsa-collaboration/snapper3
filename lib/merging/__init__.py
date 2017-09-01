"""
File contains some stuff related to merging clusters.

author: ulf.schaefer@phe.gov.uk

"""

from lib.distances import get_all_pw_dists, get_distances
from lib.ClusterMerge import ClusterMerge
from lib.ClusterStats import ClusterStats

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
        {lvl: ClusterMerge object (with lvl and org members only at this point)}
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
                oMer = ClusterMerge(level=lvl, clusters=clusters.keys())
                merges[lvl] = oMer

    return merges

# --------------------------------------------------------------------------------------------------

def get_stats_for_merge(cur, oMerge):
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
    oMerge: obj
        ClusterMerge object

    Returns
    -------
    oStats: obj
        ClusterStats object
    current_mems: list of ints
        list of members of the merged cluster
    None if problem
    """

    # get the members for each of the clusters to merge and put them in dict
    # members[clu_id] = [list of sample ids]
    clu_to_merge = oMerge.org_clusters
    members = {}
    t_lvl = oMerge.t_level
    for ctm in clu_to_merge:
        # get the members for all cluster that need merging, ignoring the clusters that fail zscore
        sql = "SELECT c.fk_sample_id FROM sample_clusters c, samples s WHERE c."+t_lvl+"=%s AND s.pk_id=c.fk_sample_id AND s.ignore_zscore IS FALSE;"
        cur.execute(sql, (ctm, ))
        rows = cur.fetchall()
        ctm_mems = [r['fk_sample_id'] for r in rows]
        members[ctm] = ctm_mems

    # this now has sample_id of the largest cluster first
    clu_to_merge = sorted(members, key=lambda k: len(members[k]), reverse=True)
    oMerge.final_name = clu_to_merge[0]

    # get stats for the biggest cluster from database
    sql = "SELECT nof_members, nof_pairwise_dists, mean_pwise_dist, stddev FROM cluster_stats WHERE cluster_level=%s AND cluster_name=%s"
    cur.execute(sql, (t_lvl, clu_to_merge[0], ))
    if cur.rowcount != 1:
        return None, None
    row = cur.fetchone()

    if row['nof_members'] > 1:

        # make a stats obj from the stuff for the biggest cluster from the db
        oMerge.stats = ClusterStats(members=row['nof_members'], stddev=row['stddev'], mean=row['mean_pwise_dist'])

        # get members of biggest cluster
        current_mems = members[clu_to_merge[0]]
        # get all the samples that need adding to it to facilitate the merge
        new_members = []
        for ctm in clu_to_merge[1:]:
            new_members += members[ctm]

        # get all distances for new members and update stats obj iteratively
        for nm in new_members:
            dists = get_distances(cur, nm, current_mems)
            all_dists_to_new_member = [d for (s, d) in dists]
            oMerge.stats.add_member(all_dists_to_new_member)
            current_mems.append(nm)

    else:
        # if the biggest cluster has only one member, get all members of all clusters to be merged
        # and get all pw distances - shouldn't be many

        # make a flat list out of the values in members which are lists
        current_mems = [item for sublist in members.values() for item in sublist]

        all_pw_dists = get_all_pw_dists(cur, current_mems)
        oMerge.stats = ClusterStats(members=len(current_mems), dists=all_pw_dists)

    oMerge.final_members = current_mems

    return oMerge.final_members

# --------------------------------------------------------------------------------------------------

def do_the_merge(cur, oMerge):
    """
    Merge the clusters on level lvl.

    Parameters
    ----------
    cur: obj
        database cursor
    oMerge: obj
        ClusterMerge object

    Returns
    -------
    None
    """

    # if we have checked z-scores we already have set the final name and calculated stats for the merge
    # if this was deactivated by the user we need to do it now.
    if oMerge.final_name == None:
        # this calculates ClusterStats for the merged cluster
        _ = get_stats_for_merge(cur, oMerge)

        # we still need to calculate the mean distance of all members of the merged cluster to all other members
        oMerge.calculate_per_member_stats(cur)

    oMerge.update_tables(cur)

    return None

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
