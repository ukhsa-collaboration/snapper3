import sys
import os
import argparse
import logging
import psycopg2
from psycopg2.extras import DictCursor
from math import sqrt

from time import time
from operator import itemgetter

__version__ = '0.1'
__date__ = '30Sep2016'
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

    parser.add_argument("--connstring",
                        "-c",
                        type=str,
                        metavar="CONNECTION",
                        required=True,
                        dest="db",
                        help="Connection string for db")

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

    logging.basicConfig(format="[%(asctime)s] %(levelname)s: %(message)s", level=logging.INFO)

    try:
        # open source db
        conn = psycopg2.connect(oArgs.db)
        cur = conn.cursor(cursor_factory=psycopg2.extras.DictCursor)

        sample_id = 1
        chr_id = 1

        t0 = time()
        cur.callproc("get_all_distances_by_id", [sample_id , chr_id])
        result = cur.fetchall()
        t1 = time()
        logging.info("Calculated %i distances with 'get_all_distances_by_id' in %.3f seconds", len(result), t1 - t0)

        t0 = time()
        cur.callproc("get_pairwise_distance", [chr_id, sample_id, 2])
        result = cur.fetchall()
        t1 = time()
        logging.info("Calculated 1 distance with 'get_pairwise_distance' in %.3f seconds", t1 - t0)

        t0 = time()
        cur.callproc("get_sample_distances_by_id", [sample_id, chr_id, range(2, 102)])
        result = cur.fetchall()
        for x in sorted(result, key=itemgetter(0)):
            print x
        t1 = time()
        logging.info("Calculated %i distances with 'get_sample_distances_by_id' in %.3f seconds", len(result), t1 - t0)

        conn.commit()

    except SystemExit as e:
        logging.error("Could not complete migration because: %s" % (str(e)))
    except psycopg2.Error as e:
         logging.error("Database reported error: %s" % (str(e)))
    finally:
        # close all dbs
        cur.close()
        conn.close()

    return 0

# end of main --------------------------------------------------------------------------------------

if __name__ == '__main__':
    sys.exit(main())

"""
chr_id = 1
samples_in = range(2, 11)

# variants = {"A":[], "C":[], "G":[], "T":[], "N":[], "-":[]}

sql = "SELECT a_pos, c_pos, g_pos, t_pos, n_pos, gap_pos FROM variants WHERE fk_sample_id=1"
cur.execute(sql)
row = cur.fetchone()
print row['a_pos']

variants = {"A": row['a_pos'], "C": row['c_pos'], "G": row['g_pos'], "T": row['t_pos'], "N": row['n_pos'], "-": row['gap_pos']}

print variants['-']
"""
