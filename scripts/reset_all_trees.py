#!/usr/bin/env python

import sys
import os
import argparse
import logging
from datetime import datetime
import psycopg2
from psycopg2.extras import DictCursor

from lib.distances import get_distances
from lib.SnapperDBInterrogation import SnapperDBInterrogation, SnapperDBInterrogationError

__version__ = '0.1'
__date__ = '17Jan2018'
__author__ = 'ulf.schaefer@phe.gov.uk'

# --------------------------------------------------------------------------------------------------

def parse_args():
    """
    Parge arguments
    Parameters
    ----------
    no inputs
    Returns
    -------
    oArgs: obj
        arguments object
    """

    sDescription = 'version %s, date %s, author %s' %(__version__, __date__, __author__)

    parser = argparse.ArgumentParser(description=sDescription)

    parser.add_argument("--debug",
                      action="store_true",
                      help="More verbose logging (default: turned off).")

    parser.add_argument("--connstring",
                        "-c",
                        type=str,
                        metavar="CONNECTION",
                        required=True,
                        dest="db",
                        help="Connection string for the db. [REQUIRED]")

    parser.add_argument("--reset",
                        "-R",
                        action="store_true",
                        default=False,
                        dest="reset",
                        help="""All trees in the db will be DELETED and recalculated. Trees table is dropped
 and recreated. Use this paramaters if you're sure about this. [REQUIRED]""")

    parser.add_argument("--reference",
                        "-r",
                        type=str,
                        required=True,
                        dest="ref",
                        help="Path to reference specified (FASTA). Required for ML, else ignored.")

    parser.add_argument("--refname",
                        type=str,
                        required=True,
                        help="The name of the reference in the database. Required for ML, else ignored.")

    oArgs = parser.parse_args()
    return oArgs

# --------------------------------------------------------------------------------------------------

def main():
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

    oArgs = parse_args()

    lglvl = logging.DEBUG if oArgs.debug else logging.INFO
    logging.basicConfig(format="[%(asctime)s] %(levelname)s: %(message)s", level=lglvl)

    if oArgs.reset == False:
        logging.error("Plese use -R flag to indicate you're sure about this.")
        return 1

    try:
        # open db
        conn = psycopg2.connect(oArgs.db)
        cur = conn.cursor(cursor_factory=psycopg2.extras.DictCursor)

        empty_trees_table(cur)

        treeables = set()

        sql = "SELECT cluster_name FROM cluster_stats WHERE nof_members>=%s AND cluster_level=%s"
        cur.execute(sql, (3, 't5', ))
        rows = cur.fetchall()
        treeables.update([row['cluster_name'] for row in rows])

        logging.info("Calculating trees for %i t5 clusters.", len(treeables))

        for t5_name in treeables:
            tree_samples, t50_size = get_tree_samples_set(cur, t5_name)

            logging.info("The tree for t5 cluster %s will contain %i samples.", t5_name, len(tree_samples))

            sql = "SELECT sample_name FROM samples WHERE pk_id IN %s"
            cur.execute(sql, (tuple(tree_samples), ))
            rows = cur.fetchall()
            sample_names = set([r['sample_name'] for r in rows])

            try:
                sample_names.remove(oArgs.refname)
            except KeyError:
                pass

            nwktree = None
            with SnapperDBInterrogation(conn_string=oArgs.db) as sdbi:
                try:
                    nwktree = sdbi.get_tree(list(sample_names),
                                           None,
                                           'ML',
                                           ref=oArgs.ref,
                                           refname=oArgs.refname,
                                           rmref=True)
                except SnapperDBInterrogationError as e:
                    logging.error(e)
                else:
                    logging.info("Tree calculation completed successfully.")

            nownow = datetime.now()

            sql = "INSERT INTO trees (nwkfile, t5_name, t50_size, sample_set, mod_date, created_at, lockdown) VALUES (%s, %s, %s, %s, %s, %s, %s)"
            cur.execute(sql, (nwktree, t5_name, t50_size, list(tree_samples), nownow, nownow, False))

        conn.commit()

    except SystemExit as e:
        logging.error("Could not complete tree creation because: %s" % (str(e)))
    except psycopg2.Error as e:
         logging.error("Database reported error: %s" % (str(e)))
    finally:
        # close all dbs
        cur.close()
        conn.close()

    return 0

# end of main --------------------------------------------------------------------------------------

def get_tree_samples_set(cur, t5_name):
    """
    Find out which samples need to go into the tree for this t5 cluster

    Parameters
    ----------
    cur: obj
        database cursor
    t5_name: int
        the name of the t5 cluster

    Returns
    -------
    sample_set: set
        it's in the name
    t50_size: int
        nof members in the t50 cluster
    """

    logging.debug("Processing t5 cluster %s", t5_name)

    sample_set = set()
    t50_cluster = set()
    t5_members = []
    t50_members = set()

    # get cluster members WITHOUT 'known outlier'
    sql = "select c.fk_sample_id AS samid, c.t50 AS tfifty FROM sample_clusters c, samples s WHERE c.t5=%s AND c.fk_sample_id=s.pk_id AND s.ignore_zscore=FALSE"
    cur.execute(sql, (t5_name, ))
    rows = cur.fetchall()
    t5_members = [r['samid'] for r in rows]
    t50_cluster.update([r['tfifty'] for r in rows])

    logging.debug("t5 %s has %i members.", t5_name, len(t5_members))

    # all members of the t5 are definitely in the tree
    sample_set.update(t5_members)

    assert len(t50_cluster) == 1, ("Why are not all members of t5 %s in the same t50?" % (t5_name))

    t50_name = t50_cluster.pop()

    logging.debug("t5 %s sits within t50 %s", t5_name, t50_name)

    sql = "select c.fk_sample_id AS samid FROM sample_clusters c, samples s WHERE c.t50=%s AND c.fk_sample_id=s.pk_id AND s.ignore_zscore=FALSE"
    cur.execute(sql, (t50_name, ))
    rows = cur.fetchall()
    t50_members.update([r['samid'] for r in rows])
    t50_size = len(t50_members)

    logging.debug("t50 %s has %i members.", t50_name, t50_size)

    for t5_mem in t5_members:
        check_samples = t50_members.difference(sample_set)
        dists = get_distances(cur, t5_mem, list(check_samples))
        sample_set.update([sa for (sa, di) in dists if di <= 50])

    return sample_set, t50_size

# -------------------------------------------------------------------------------------------------

def empty_trees_table(cur):
    """
    Drop and recreate the trees table

    Parameters
    ----------
    cur: obj
        database cursor

    Returns
    -------
    always 0
    """

    sql = "DROP TABLE trees"
    cur.execute(sql)

    sql = "CREATE TABLE trees (pk_id SERIAL PRIMARY KEY, nwkfile bytea, t5_name integer, t50_size integer, sample_set integer[], mod_date timestamp, created_at timestamp, lockdown boolean DEFAULT FALSE)"
    cur.execute(sql)

    return 0

# -------------------------------------------------------------------------------------------------

if __name__ == '__main__':
    sys.exit(main())
