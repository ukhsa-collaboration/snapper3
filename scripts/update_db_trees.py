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
            sql = "SELECT pk_id, sample_set, mod_date, t50_size FROM trees WHERE t5_name=%s"
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
                mod_time = row['mod_date']
                t50_size = row['t50_size']

                update_an_existing_tree(cur,
                                        conn,
                                        tree_row_id,
                                        t5_c,
                                        t5_mems,
                                        tree_sample_set,
                                        mod_time,
                                        t50_size,
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
        sample_names.remove(args['refname'])
    except KeyError:
        pass

    nwktree = None
    with SnapperDBInterrogation(conn_string=args['db']) as sdbi:
        try:
            nwktree = sdbi.get_tree(list(sample_names),
                                   None,
                                   'ML',
                                   ref=args['ref'],
                                   refname=args['refname'],
                                   rmref=True)
        except SnapperDBInterrogationError as e:
            logging.error(e)
        else:
            logging.info("Tree calculation completed successfully.")

    nownow = datetime.now()

    sql = "INSERT INTO trees (nwkfile, t5_name, t50_size, sample_set, mod_date, created_at, lockdown) VALUES (%s, %s, %s, %s, %s, %s, %s)"
    cur.execute(sql, (nwktree, t5_name, len(t50_members), list(sample_set), nownow, nownow, False, ))

    return 0

# --------------------------------------------------------------------------------------------------

def update_an_existing_tree(cur, conn, tree_row_id, t5_name, t5_members, tree_sample_set, mod_time, t50_size, args):
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
    mod_time: datetime.datetime
        time tree was last updated
    t50_size: int
        size of the t50 cluster last time it was updated
    args: dict
        as passed to main function

    Returns
    -------
    1 or None if fail
    """

    logging.debug("=== Checking if tree for t5 cluster %s needs updating. ===", t5_name)

    t50_name = get_t50_cluster(cur, t5_name, t5_members)
    logging.debug("t5 %s sits within t50 %s. It has %i members", t5_name, t50_name, len(t5_members))

    t50_members = get_members(cur, 't50', t50_name)
    logging.debug("t50 %s has %i members.", t50_name, len(t50_members))

    logging.debug("t50 size at last update was %i, t50 size now is %i", t50_size, len(t50_members))
    if len(t50_members) <= t50_size:
        logging.debug("Tree for t5 cluster %s does not need updating.", t5_name)
        return 1

    needs_update = False

    # get the maximum t0 cluster number in the previous tree
    sql = "SELECT max(t0) FROM sample_clusters WHERE fk_sample_id IN %s"
    cur.execute(sql, (tuple(tree_sample_set), ), )
    tree_t0_max = cur.fetchone()[0]

    logging.debug("Max t0 in this tree is: %i", tree_t0_max)

    # set with t5 members that are not in the tree yet
    new_t5_members = t5_members.difference(tree_sample_set)
    # set with t5 members already in the tree
    old_t5_members = t5_members.intersection(tree_sample_set)
    # for all t5 members that are not in the tree yet
    # check whether there are members in the t50 that are not in the tree yet, but should be
    logging.debug("There are %i new members in this t5 cluster.", len(new_t5_members))
    for new_t5_mem in new_t5_members:
        # definitely needs updating when there is new t5 members
        needs_update = True
        # check only the distances to t50 cluster samples that are not in the tree yet
        check_samples = t50_members.difference(tree_sample_set)
        logging.debug("Checking distances from new t5 member %s to %i t50 members that are not in the tree yet.", new_t5_mem, len(check_samples))
        dists = get_distances(cur, new_t5_mem, list(check_samples))
        new_members = [sa for (sa, di) in dists if di <= 50]
        if len(new_members) > 0:
            logging.debug("These samples need to be in the tree now: %s.", str(new_members))
            tree_sample_set.update(new_members)

    # reduce the list of members in the 50 to those REALLY neeeding checking
    filtered_t50_members = filter_samples_to_be_checked(cur, t50_members, tree_t0_max)
    logging.debug("Reduced nof t50 members to be checked from %i to %i.", len(t50_members), len(filtered_t50_members))

    for old_t5_mem in old_t5_members:
        # check only the distances to t50 cluster samples that are not in the tree yet
        check_samples = filtered_t50_members.difference(tree_sample_set)
        logging.debug("Checking distances from old t5 member %s to %i t50 members that are not in the tree yet.", old_t5_mem, len(check_samples))
        dists = get_distances(cur, old_t5_mem, list(check_samples))
        new_members = [sa for (sa, di) in dists if di <= 50]
        if len(new_members) > 0:
            logging.debug("These samples need to be in the tree now: %s.", str(new_members))
            tree_sample_set.update(new_members)
            needs_update = True

    if needs_update == True:

        # lock table row during tree update
        sql = "UPDATE trees SET nwkfile=%s, sample_set=%s, lockdown=%s WHERE pk_id=%s"
        cur.execute(sql, (None, None, True, tree_row_id, ))
        conn.commit()

        logging.info("The tree for t5 cluster %s needs updating and will now contain %i samples.", t5_name, len(tree_sample_set))

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
                                        refname=args['refname'],
                                        rmref=True)
            except SnapperDBInterrogationError as e:
                logging.error(e)
                return None
            else:
                logging.info("Tree calculation completed successfully.")

        nownow = datetime.now()

        # update the database - unlock the table row
        sql = "UPDATE trees SET nwkfile=%s, lockdown=%s, mod_date=%s, sample_set=%s, t50_size=%s WHERE pk_id=%s"
        cur.execute(sql, (nwktree, False, nownow, list(tree_sample_set), len(t50_members), tree_row_id, ))
        conn.commit()

    else:
        sql = "UPDATE trees SET t50_size=%s WHERE pk_id=%s"
        cur.execute(sql, (len(t50_members), tree_row_id, ))
        logging.debug("Tree for t5 cluster %i does not need updating.", t5_name)

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

def filter_samples_to_be_checked(cur, samples, tree_t0_max):
    """
    Get the name of the t50 cluster that all members of this t5 cluster a members of.

    Parameters
    ----------
    cur: obj
        database cursor object
    samples: set
        unfiltered set of samples
    tree_t0_max: int
        max t0 of any previous tree members

    Returns
    -------
    x: set
        minimal set of samples that need to be checked
    """

    t0s = {}
    t0len = {}

    # get the t0 for each member in the set of samples
    sql = "SELECT fk_sample_id, t0 FROM sample_clusters WHERE fk_sample_id IN %s"
    cur.execute(sql, (tuple(samples), ))
    rows = cur.fetchall()

    for r in rows:
        t0s[r['fk_sample_id']] = r['t0']
        try:
            _ = t0len[r['t0']]
        except KeyError:
            sql = "SELECT fk_sample_id FROM sample_clusters WHERE t0=%s"
            cur.execute(sql, (r['t0'], ))
            cur.fetchall()
            t0len[r['t0']] = cur.rowcount

    # t0s[sam1] = 123
    # t0len[123] = 2 <- length

    # return only those samples that either:
    #     - have a t0 that is the larger than the biggest in the tree so far
    #     - are not alone in their t0
    # but why?
    # => because samples that have a smaller t0 than one that exists in the tree have already
    #    been checked against current tree samples
    #    Except: the ones that have been added to an existing t0 or are in a t0 that has been merged
    #            -> therefore filter out only t0 singletons

    ret = set([s for s in samples if (t0s[s] >= tree_t0_max or t0len[t0s[s]] >= 2)])

    logging.debug("These samples were excluded from consideration: %s", str(samples.difference(ret)))

    return ret

# --------------------------------------------------------------------------------------------------

if __name__ == "__main__":
    exit(main(vars(get_args().parse_args())))
