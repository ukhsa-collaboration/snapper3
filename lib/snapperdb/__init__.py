"""
File contains some stuff that is used in snapper v3.

author: ulf.schaefer@phe.gov.uk

"""

import logging
from time import time
from operator import itemgetter

from lib.utils import get_closest_threshold

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

def get_relevant_distances(cur, sample_id):
    """
    Get the distances to this sample from the database.

    Parameters
    ----------
    cur: obj
        database cursor
    sample_id: int
        sample pk_id

    Returns
    -------
    d: list of tuples
        sorted list of tuples with (sample_id, distance) with closes sample first
        e.g. [(298, 0), (37, 3), (55, 4)]
        None if fail
    """

    # get the relevant samples from the database, these are the ones that have been clustered and are not ignored
    sql = "SELECT c.fk_sample_id FROM sample_clusters c, samples s WHERE s.pk_id=c.fk_sample_id AND s.ignore_sample IS FALSE"
    cur.execute(sql)
    rows = cur.fetchall()
    relv_samples = [r['fk_sample_id'] for r in rows]

    # get list of contig ids from database
    sql = "SELECT pk_id FROM contigs"
    cur.execute(sql)
    rows = cur.fetchall()
    contig_ids = [r['pk_id'] for r in rows]

    d = {}
    for cid in contig_ids:
        t0 = time()
        cur.callproc("get_sample_distances_by_id", [sample_id, cid, relv_samples])
        result = cur.fetchall()
        t1 = time()
        logging.info("Calculated %i distances on contig %i with 'get_sample_distances_by_id' in %.3f seconds", len(result), cid, t1 - t0)

        # sum up if there are more than one contigs
        for res in result:
            if res[2] == None:
                res[2] = 0
            try:
                d[res[0]] += res[2]
            except KeyError:
                d[res[0]] = res[2]

    # d = sorted(d.items(), key=itemgetter(1), reverse=False)
    d = sorted([list(x) for x in d.items()], key=itemgetter(1), reverse=False)

    return d

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

# --------------------------------------------------------------------------------------------------

# --------------------------------------------------------------------------------------------------

# --------------------------------------------------------------------------------------------------

# --------------------------------------------------------------------------------------------------

# --------------------------------------------------------------------------------------------------

# --------------------------------------------------------------------------------------------------

# --------------------------------------------------------------------------------------------------

# --------------------------------------------------------------------------------------------------
