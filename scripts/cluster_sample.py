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

    logging.debug("Args received: %s", str(args))

    try:
        # open db
        conn = psycopg2.connect(args['db'])
        cur = conn.cursor(cursor_factory=psycopg2.extras.DictCursor)

        sample_id = sndb.get_sample_id(cur, args['sample_name'])
        if sample_id < 0:
            logging.error("Could not get sample id from db. :-(")
            return 1

        logging.info("Processing sample %s with id %i", args['sample_name'], sample_id)

        distances = sndb.get_relevant_distances(cur, sample_id)
        if distances == None:
            logging.error("Could not get distances from db. :-(")
            return 1







        # distances = distances[2:]
        # distances[5][1] = 5
        # distances[6][1] = 5
        #distances = [(x[0], x[1]+16)for x in distances]







        logging.debug("Distances calculated: %s", str(distances))

        nbhood = sndb.get_closest_samples(cur, distances)
        """
        nbhood = {'closest_distance': int,
                  'nearest_t': int,
                  'closest_sample': int,
                  'closest_snad': [list of 7 ints]}
        """
        logging.debug("Sample neighbourhood: %s", str(nbhood))

        new_snad = sndb.get_new_snp_address(nbhood)

        logging.debug("Proposed SNP address for this sample: %s-%s-%s-%s-%s-%s-%s",
                      new_snad[6], new_snad[5], new_snad[4], new_snad[3], new_snad[2], new_snad[1], new_snad[0])

        merges = sndb.check_merging_needed(cur, distances, new_snad)

        logging.debug("Merges that would be required to make this assignment: %s", str(merges))

        zscore_fail, zscore_info = sndb.check_zscores(cur, distances, new_snad, nbhood, merges)

        if zscore_fail == None:
            logging.error("Could not calculate z-scores. :-(")
            return 1

        if zscore_fail == True:
            logging.error("z-score check for this assignment has failed. Database is not updated.")
            for s in zscore_info:
                logging.info(s)
            return 1

        logging.debug("All z-score checks passed for this assignment.")

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
