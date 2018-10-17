"""
File contains some to get some distance calculations done for snapperdb v3.

author: ulf.schaefer@phe.gov.uk

"""

import logging
from time import time
import json
import gzip
from operator import itemgetter

import requests

# --------------------------------------------------------------------------------------------------

def get_relevant_samples(cur):
    """
    Get the relevant samples from the database, these are the ones that have been clustered and are not ignored

    Parameters
    ----------
    cur: obj
        database cursor

    Returns
    -------
    relv_samples: list
        list of ints which are the sample ids
    """
    # get the relevant samples from the database, these are the ones that have been clustered and are not ignored
    sql = "SELECT c.fk_sample_id FROM sample_clusters c, samples s WHERE s.pk_id=c.fk_sample_id AND s.ignore_sample IS FALSE"
    cur.execute(sql)
    rows = cur.fetchall()
    relv_samples = [r['fk_sample_id'] for r in rows]
    logging.info("Found %i samples in the database relavant for clustering.", len(relv_samples))
    return relv_samples

# --------------------------------------------------------------------------------------------------

def get_all_pw_dists(cur, samids):
    """
    Get all pairwise distances between the samples in the input list.

    Parameters
    ----------
    cur: obj
        database cursor
    samids: list of int
        sample ids

    Returns
    -------
    dists: lists of ints
        lists with distances
    None if there is a problem
    """

    # get list of contig ids from database
    sql = "SELECT pk_id FROM contigs"
    cur.execute(sql)
    rows = cur.fetchall()
    contig_ids = [r['pk_id'] for r in rows]

    samids = set(samids)
    done = set()
    dists = []

    for s in samids:
        # note the ones that are already done, so nothing is calculated twice.
        done.add(s)
        oths = samids.difference(done)

        d = {}
        for cid in contig_ids:

            cur.callproc("get_sample_distances_by_id", [s, cid, list(oths)])
            result = cur.fetchall()

            for res in result:
                if res[2] == None:
                    res[2] = 0
                try:
                    d[res[0]] += res[2]
                except KeyError:
                    d[res[0]] = res[2]

        dists += d.values()

    assert len(dists) == (len(samids) * (len(samids)-1))/2

    return dists

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

    relv_samples = get_relevant_samples(cur)
    d = get_distances(cur, sample_id, relv_samples)

    return d

# --------------------------------------------------------------------------------------------------

def get_missing_distances(cur, sample_id, haves):
    """
    Get the distances to this sample that we don't already have from the database.

    Parameters
    ----------
    cur: obj
        database cursor
    sample_id: int
        sample pk_id
    haves: set
        set of sample ids tuples for which we already have the distance.

    Returns
    -------
    d: list of tuples
        sorted list of tuples with (sample_id, distance) with closes sample first
        e.g. [(298, 0), (37, 3), (55, 4)]
        None if fail
    """

    all_samples = set(get_relevant_samples(cur))
    relv_samples = list(all_samples.difference(haves))
    d = get_distances(cur, sample_id, relv_samples)

    return d


# --------------------------------------------------------------------------------------------------

def get_distances(cur, samid, others):
    """
    Get the distances of this sample to the other samples from the database.

    Parameters
    ----------
    cur: obj
        database cursor
    sample_id: int
        sample pk_id
    others: list of int
        other samples to calculate the distance to

    Returns
    -------
    d: list of tuples
        sorted list of tuples with (sample_id, distance) with closes sample first
        e.g. [(298, 0), (37, 3), (55, 4)]
        None if fail
    """

    # get list of contig ids from database
    sql = "SELECT pk_id FROM contigs"
    cur.execute(sql)
    rows = cur.fetchall()
    contig_ids = [r['pk_id'] for r in rows]

    d = {}
    for cid in contig_ids:
        t0 = time()
        cur.callproc("get_sample_distances_by_id", [samid, cid, others])
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

    d = sorted(d.items(), key=itemgetter(1), reverse=False)

    return d

# --------------------------------------------------------------------------------------------------

def get_distance_matrix(cur, samids):
    """
    Get a distance matrix for the given samples.

    Parameters
    ----------
    cur: obj
        database cursor
    samids: int
        list of sample pk_ids

    Returns
    -------
    dist: dict
        complete matrix
        dist[s1][s2] = d
        dist[s2][s1] = d
    """

    # get list of contig ids from database
    sql = "SELECT pk_id FROM contigs"
    cur.execute(sql)
    rows = cur.fetchall()
    contig_ids = [r['pk_id'] for r in rows]

    samids = set(samids)
    done = set()
    dists = {}

    for s in samids:

        try:
            dists[s][s] = 0
        except KeyError:
            dists[s] = {s: 0}

        # note the ones that are already done, so nothing is calculated twice.
        done.add(s)
        oths = samids.difference(done)

        d = {}
        for cid in contig_ids:

            cur.callproc("get_sample_distances_by_id", [s, cid, list(oths)])
            result = cur.fetchall()

            for res in result:
                if res[2] == None:
                    res[2] = 0
                try:
                    d[res[0]] += res[2]
                except KeyError:
                    d[res[0]] = res[2]
        for osam, snpdi in d.items():
            dists[s][osam] = snpdi
            try:
                dists[osam][s] = snpdi
            except KeyError:
                dists[osam] = {s: snpdi}

    return dists

# --------------------------------------------------------------------------------------------------

def get_distances_precalc(cur, sam_id, sample_name, json_file_name):
    """
    Check the precalculated data against the database, get the missing distances,
    and add the precalculated ones.

    Parameters
    ----------
    cur: obj
        database cursor
    sample_id: int
        sample pk_id
    sample_name: str
        samples name
    json_file_name: str
        name of a json file containing

    Returns
    -------
    d: list of tuples
        sorted list of tuples with (sample_id, distance) with closes sample first
        e.g. [(298, 0), (37, 3), (55, 4)]
        None if fail
    """

    open_func = gzip.open if json_file_name.endswith('.gz') == True else open
    precalc_data = None
    with open_func(json_file_name) as f:
        precalc_data = json.load(f)

    # chck the sample name in the json against the current sample
    if precalc_data['sample_name'] != sample_name:
        logging.error("Sample name does not match precalculated data!")
        return None

    # check data consistency between the db and the precalculated distances, do samid and samname match?
    sql = "SELECT c.fk_sample_id AS sid, s.sample_name AS name FROM sample_clusters c, samples s WHERE s.pk_id=c.fk_sample_id AND s.ignore_sample IS FALSE"
    cur.execute(sql)
    rows = cur.fetchall()
    tbl_values = {r['sid']: r['name'] for r in rows}
    for (pre_id, pre_name, dis) in precalc_data['distances']:
        if tbl_values[pre_id] != pre_name:
            logging.error("Precalculated data does not match database. Precalculated samples name for id %i was %s, but in db it's %s",
                          pre_id, pre_name, tbl_values[pre_id])
            return None

    # get missing distances, add the precalculated ones and re-sort
    d = get_missing_distances(cur, sam_id, set([x[0] for x in precalc_data['distances']]))
    d += [(x[0], x[2]) for x in precalc_data['distances']]
    d.sort(key=lambda x: x[1])

    return d

# --------------------------------------------------------------------------------------------------

def get_distances_fusion(cur, sample_id, fusion_url):
    """
    Get the relevant distances back from the fusion web service.

    Parameters
    ----------
    cur: obj
        database cursor
    sample_id: int
        sample pk_id
    fusion_url: str
        the url of the fusion webservice master

    Returns
    -------
    d: list of tuples
        sorted list of tuples with (sample_id, distance) with closes sample first
        e.g. [(298, 0), (37, 3), (55, 4)]
        None if fail
    """

    relv_samples = get_relevant_samples(cur)

    tic = time()

    res = requests.post('%s/some_distances/%s' % (fusion_url, sample_id),
                        json=[str(x) for x in relv_samples],
                        headers={'Cache-Control': 'no-cache'})

    toc = time()

    if res.status_code != 200:
        logging.error("There was a problem getting the distances from fusion.")
        logging.error("Please consult fusion server logs.")
        return None

    data = res.json()
    logging.info("%i distances for sample %s were successfully calculated by fusion webservice %s in %.3f secs.",
                 len(data['distances']),
                 sample_id,
                 fusion_url,
                 toc - tic)

    if len(data['missing_samples']) > 0:
        logging.error("The distances to the following samples could not be calculated: %s",
                      str(data['missing_samples']))

    return [(int(x[0]), x[1]) for x in data['distances']]

# --------------------------------------------------------------------------------------------------
