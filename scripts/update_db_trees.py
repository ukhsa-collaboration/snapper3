import logging
import sys
import os
from datetime import datetime
import argparse
from argparse import RawTextHelpFormatter

import psycopg2
from psycopg2.extras import DictCursor

from lib.distances import get_distances
from lib.SnapperDBInterrogation import SnapperDBInterrogation, SnapperDBInterrogationError

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

    return r'''Check if there are any samples in the database that are not in the t5_tree they should be in (if any)
and recalculate those trees.'''

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

    args.add_argument("--reference",
                        "-r",
                        type=str,
                        required=True,
                        dest="ref",
                        help="Path to reference specified REQUIRED.")

    args.add_argument("--refname",
                        type=str,
                        required=True,
                        help="The name of the reference in the database. REQUIRED.")

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

        # get all t5 clusters that have 3 or more members
        t5_clusters = list()
        sql = "SELECT cluster_name, nof_members FROM cluster_stats WHERE cluster_level=%s"
        cur.execute(sql, ('t5', ))
        rows = cur.fetchall()
        for r in rows:
            t5_clusters = [r['cluster_name'] for r in rows if r['nof_members'] >= 3]

        logging.info("Checking trees for %i t5 clusters with >=3 members.", len(t5_clusters))

        for t5_c in t5_clusters:
            # get the members of this t5 cluster
            t5_mems = get_members(cur, 't5', t5_c)

            # is there a tree for this?
            sql = "SELECT pk_id, sample_set FROM trees WHERE t5_name=%s"
            cur.execute(sql, (t5_c, ))
            if cur.rowcount == 0:
                # no
                logging.info("Making a new tree for t5 cluster %i", t5_c)
                make_a_new_tree(cur, t5_c, t5_mems, args)

            elif cur.rowcount == 1:
                # yes
                row = cur.fetchone()

                # get the pkid and the set of sample ids in this tree
                tree_row_id = row['pk_id']
                tree_sample_set = set(row['sample_set'])

                # are all samples that are currently in the cluster (including recent additions)
                # already in the tree?
                if t5_mems.issubset(tree_sample_set) == True:
                    # yes
                    logging.info("Tree for t5 cluster %i is still up to date.", t5_c)
                else:
                    # no
                    logging.info("Updating an existing tree for t5 cluster %i in rows with pk_id %s", t5_c, tree_row_id)
                    logging.info("This is because samples %s should be in the tree but are not.", t5_mems.difference(tree_sample_set))
                    update_an_existing_tree(cur,
                                            conn,
                                            tree_row_id,
                                            t5_c,
                                            t5_mems,
                                            tree_sample_set,
                                            args)

            else:
                logging.error("Multiple trees for t5 cluster %i", t5_c)

        conn.commit()

    except psycopg2.Error as e:
         logging.error("Database reported error: %s" % (str(e)))
    finally:
        # close all dbs
        cur.close()
        conn.close()

    return 0

# end of main --------------------------------------------------------------------------------------

def make_a_new_tree(cur, t5_name, t5_members, args):
    """
    Makes a new tree for a cluster that previously did not have one.

    Parameters
    ----------
    t5_name: int
        name of the t5 cluster
    t5_members: set
        set of members of the t5 cluster
    args: dict
        as passed to main function

    Returns
    -------
    0
    """

    t50_name = get_t50_cluster(cur, t5_name, t5_members)
    logging.debug("t5 %s sits within t50 %s", t5_name, t50_name)

    t50_members = get_members(cur, 't50', t50_name)
    logging.debug("t50 %s has %i members.", t50_name, len(t50_members))

    sample_set = set()
    sample_set.update(list(t5_members))

    for t5_mem in t5_members:
        check_samples = t50_members.difference(sample_set)
        dists = get_distances(cur, t5_mem, list(check_samples))
        sample_set.update([sa for (sa, di) in dists if di <= 50])

    logging.info("The tree for t5 cluster %s will contain %i samples.", t5_name, len(sample_set))

    sample_names = get_sample_names(cur, sample_set)

    try:
        sample_names.remove(oArgs.refname)
    except KeyError:
        pass

    nwktree = None
    with SnapperDBInterrogation(conn_string=args['db']) as sdbi:
        try:
            nwktree = sdbi.get_tree(list(sample_names),
                                   None,
                                   'ML',
                                   ref=args['ref'],
                                   refname=args['refname'])
        except SnapperDBInterrogationError as e:
            logging.error(e)
        else:
            logging.info("Tree calculation completed successfully.")

    nownow = datetime.now()

    sql = "INSERT INTO trees (nwkfile, t5_name, sample_set, mod_date, created_at, lockdown) VALUES (%s, %s, %s, %s, %s, %s)"
    cur.execute(sql, (nwktree, t5_name, list(sample_set), nownow, nownow, False))

    return 0

# --------------------------------------------------------------------------------------------------

def update_an_existing_tree(cur, conn, tree_row_id, t5_name, t5_members, tree_sample_set, args):
    """
    Updates an existing tree in the database

    Parameters
    ----------
    cur: obj
        database cursor object
    conn: onj
        database connection object
    tree_row_id: int
        pk_id in the trees table
    t5_name: int
        name of the t5 cluster
    t5_members: set
        set of members for this t5 cluster
    tree_sample_set: set
        set of samples in the tree
    args: dict
        as passed to main function

    Returns
    -------
    1 or None if fail
    """

    sql = "UPDATE trees SET nwkfile=%s, sample_set=%s, lockdown=%s WHERE pk_id=%s"
    cur.execute(sql, (None, None, True, tree_row_id, ))
    conn.commit()

    tree_sample_set.update(t5_members)
    logging.debug("t5 %s has %i members.", t5_name, len(t5_members))

    t50_name = get_t50_cluster(cur, t5_name, t5_members)
    logging.debug("t5 %s sits within t50 %s", t5_name, t50_name)

    t50_members = get_members(cur, 't50', t50_name)
    logging.debug("t50 %s has %i members.", t50_name, len(t50_members))

    for t5_mem in t5_members:
        # check only the distances to t50 cluster samples that are not in the tree yet
        check_samples = t50_members.difference(tree_sample_set)
        dists = get_distances(cur, t5_mem, list(check_samples))
        tree_sample_set.update([sa for (sa, di) in dists if di <= 50])

    logging.info("The tree for t5 cluster %s will contain %i samples.", t5_name, len(tree_sample_set))

    sample_names = get_sample_names(cur, tree_sample_set)

    # if the reference is part of the tree we need to remove this here
    # it is always part of all trees anyway
    try:
        sample_names.remove(args['refname'])
    except KeyError:
        pass

    # make a tree now using SnapperDBInterrogation interface
    nwktree = None
    with SnapperDBInterrogation(conn_string=args['db']) as sdbi:
        try:
            nwktree = sdbi.get_tree(list(sample_names),
                                   None,
                                   'ML',
                                   ref=args['ref'],
                                   refname=args['refname'])
        except SnapperDBInterrogationError as e:
            logging.error(e)
            return None
        else:
            logging.info("Tree calculation completed successfully.")

    nownow = datetime.now()

    # update the databse
    sql = "UPDATE trees SET nwkfile=%s, lockdown=%s, mod_date=%s, sample_set=%s WHERE pk_id=%s"
    cur.execute(sql, (nwktree, False, nownow, list(tree_sample_set), tree_row_id, ))
    conn.commit()

    return 1

# --------------------------------------------------------------------------------------------------

def get_sample_names(cur, tree_sample_set):
    """
    Get the set of members for a given t5 cluster, without 'known outliers'

    Parameters
    ----------
    cur: obj
        database cursor object
    tree_sample_set: set or list
        interable with all the sample ids in the tree

    Returns
    -------
    sample_names: set
        names of these samples
    """

    sql = "SELECT sample_name FROM samples WHERE pk_id IN %s"
    cur.execute(sql, (tuple(tree_sample_set), ))
    rows = cur.fetchall()
    sample_names = set([r['sample_name'] for r in rows])

    return sample_names

# --------------------------------------------------------------------------------------------------

def get_members(cur, level, name):
    """
    Get the set of members for a given t5 cluster, without 'known outliers'

    Parameters
    ----------
    cur: obj
        database cursor object
    level: str
        e.g. 't5'
    name: str
        t5 cluster name

    Returns
    -------
    t5_mems: set
        members
    """

    mems = set()
    sql  = "SELECT c.fk_sample_id FROM sample_clusters c, samples s WHERE "+level+"=%s AND c.fk_sample_id=s.pk_id AND s.ignore_zscore=FALSE"
    cur.execute(sql, (name, ))
    rows = cur.fetchall()
    mems.update([r['fk_sample_id'] for r in rows])

    return mems

# --------------------------------------------------------------------------------------------------

def get_t50_cluster(cur, t5_name, t5_members):
    """
    Get the name of the t50 cluster that all members of this t5 cluster a members of.

    Parameters
    ----------
    cur: obj
        database cursor object
    t5_name: int
        name of t5 cluster
    t5_members: list or set
        an interable with all t5 member names

    Returns
    -------
    t50_name: int
        the name of the 50 cluster
    """

    t50_cluster = set()

    sql = "SELECT c.t50 AS tfifty FROM sample_clusters c, samples s WHERE c.fk_sample_id IN %s AND c.fk_sample_id=s.pk_id AND s.ignore_zscore=FALSE"
    cur.execute(sql, (tuple(t5_members), ))
    rows = cur.fetchall()
    t50_cluster.update([r['tfifty'] for r in rows])

    assert len(t50_cluster) == 1, ("Why are not all members of t5 %s in the same t50?" % (t5_name))

    t50_name = t50_cluster.pop()

    return t50_name

# --------------------------------------------------------------------------------------------------

if __name__ == "__main__":
    exit(main(vars(get_args().parse_args())))
