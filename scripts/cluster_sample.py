import logging

import argparse
from argparse import RawTextHelpFormatter
import psycopg2
from psycopg2.extras import DictCursor

import lib.snapperdb as sndb

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

    return r'''After the variants for a sample have been added to the database, use this to
determine the clustering for this sample. Will perform all statictical checks and
merging if necessary and update the database accordingly. If statistical checks fail
database will not be updated.'''

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

    args.add_argument("--sample-name",
                      "-s",
                      type=str,
                      metavar="NAME",
                      required=True,
                      dest="sample_name",
                      help="The name of the sample to go into the db. REQUIRED.")

    args.add_argument("--justgetonwithit",
                      action='store_true',
                      help="Do not perform checks and just add the sample. It's fine. [Default: Perform checks.]")


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

    print args

    try:
        # open db
        conn = psycopg2.connect(args['db'])
        cur = conn.cursor(cursor_factory=psycopg2.extras.DictCursor)

        sample_id = sndb.get_sample_id(cur, args['sample_name'])
        if sample_id < 0:
            logging.error("Could not get sample id from db")
            return 1

        logging.info("Processing sample %s with id %i", args['sample_name'], sample_id)

        distances = sndb.get_relevant_distances(cur, sample_id)
        if distances == None:
            logging.error("Could not get distances from db")
            return 1


        conn.commit()

    except psycopg2.Error as e:
         logging.error("Database reported error: %s" % (str(e)))
    finally:
        # close all dbs
        cur.close()
        conn.close()

    return 0

# --------------------------------------------------------------------------------------------------

if __name__ == "__main__":
    exit(main(vars(get_args().parse_args())))
