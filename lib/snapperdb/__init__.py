"""
File contains some stuff that is used in snapper v3.

author: ulf.schaefer@phe.gov.uk

"""

import logging
from time import time
from operator import itemgetter

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

    d = sorted(d.items(), key=itemgetter(1), reverse=False)

    return d

# --------------------------------------------------------------------------------------------------
