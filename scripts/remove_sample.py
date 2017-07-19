import logging

import argparse
from argparse import RawTextHelpFormatter
import psycopg2
from psycopg2.extras import DictCursor

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

    return r'''This removes a samples from the database, which is not a trivial thing to do because
cluster stats and stats for other samples need to be updated.'''

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

    args.add_argument("--sample",
                      "-s",
                      type=str,
                      required=True,
                      dest="sample",
                      help="Name of sample to remove. REQUIRED.")

    args.add_argument("--just-ignore",
                      action='store_true',
                      help="""Sample and variant information will be retained in database,
but clustering information will be removed. Ignore_sample will
be set to TRUE. [DEFAULT: Remove everything. Sample can be added
and clustered again later.]""")

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

        sql = "SELECT pk_id, ignore_sample FROM samples WHERE sample_name=%s"
        cur.execute(sql, (args['sample'], ))
        if cur.rowcount < 1:
            logging.error("Sample '%s' not found in database.", args['sample'])
            return 1
        row = cur.fetchone()

        sample_id = row['pk_id']
        igno_flag = row['ignore_sample']

        if igno_flag == True:
            logging.info("This sample is already ignored.")
            if args['just_ignore'] == True:
                logging.info("You chose not to remove it completely.")
                return 0
            else:
                drop_sample(cur, sample_id)
                conn.commit()
                return 0


        sql = "SELECT t0, t5, t10, t25, t50, t100, t250 FROM sample_clusters WHERE fk_sample_id=%s"
        cur.execute(sql, (sample_id, ))
        if cur.rowcount != 1:
            logging.error("There is not exactly one entry in sample clusters for this samples. :-(")
            return 1
        row = cur.fetchone()
        snad = [row['t0'], row['t5'], row['t10'], row['t25'], row['t50'], row['t100'], row['t250']]

        print snad


        conn.commit()

    except psycopg2.Error as e:
         logging.error("Database reported error: %s" % (str(e)))
         return None
    finally:
        # close all dbs
        cur.close()
        conn.close()

    return 0

# --------------------------------------------------------------------------------------------------

def drop_sample(cur, sid):
    """
    Remove the sample with the id from the variants and samples tables.

    Parameters
    ----------
    cur: obj
        database cursor
    sid: int
        pk id of sample to remove

    Returns
    -------
    0
    """

    logging.info("Removing the sample from the samples and variants table.")

    sql = "DELETE FROM variants WHERE fk_sample_id=%s"
    cur.execute(sql, (sid, ))

    sql = "DELETE FROM samples WHERE pk_id=%s"
    cur.execute(sql, (sid, ))

    return 0

# --------------------------------------------------------------------------------------------------

if __name__ == "__main__":
    exit(main(vars(get_args().parse_args())))
