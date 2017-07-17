import logging

import argparse
from argparse import RawTextHelpFormatter
import psycopg2
from psycopg2.extras import DictCursor

import lib.snapperdb as sndb
import lib.registration as regis
import lib.merging as merging
from lib.distances import get_all_pw_dists, get_relevant_distances

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

    args.add_argument("--no-zscore-check",
                      action='store_true',
                      help="""Do not perform checks and just add the sample. It's fine.
[Default: Perform checks.]""")

    args.add_argument("--with-registration",
                      action='store_true',
                      help="""Register the clustering for this sample in the database
and update the cluster stats. [Default: Do not register.]""")

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

        distances = get_relevant_distances(cur, sample_id)
        if distances == None:
            logging.error("Could not get distances from db. :-(")
            return 1

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

        logging.info("Proposed SNP address for this sample: %s-%s-%s-%s-%s-%s-%s",
                      new_snad[6], new_snad[5], new_snad[4], new_snad[3], new_snad[2], new_snad[1], new_snad[0])

        merges = merging.check_merging_needed(cur, distances, new_snad)

        logging.info("Merges that would be required to make this assignment: %s", str([str(m) for m in merges.values()]))

        if args['no_zscore_check'] == False:
            zscore_fail, zscore_info = sndb.check_zscores(cur, distances, new_snad, nbhood, merges)

            if zscore_fail == None:
                logging.error("Could not calculate z-scores. :-(")
                return 1

            if zscore_fail == True:
                logging.error("z-score check for this assignment has failed. Database is not updated.")
                for s in zscore_info:
                    logging.info(s)
                return 1

            logging.info("All z-score checks passed for this assignment.")
        else:
            logging.info("User disabled z-score checks for this assignment.")

        logging.debug("Merges that would be required to make this assignment: %s", str([str(m) for m in merges.values()]))

        if args['with_registration'] == True:

            for lvl in merges.keys():
                merging.do_the_merge(cur, merges[lvl])

            final_snad = regis.register_sample(cur, sample_id, distances, new_snad)

            if final_snad != None:
                logging.info("Sample with sample_id %s was registered in the database with SNP address: %s-%s-%s-%s-%s-%s-%s",
                             sample_id, final_snad[6], final_snad[5], final_snad[4], final_snad[3], final_snad[2], final_snad[1], final_snad[0])
            else:
                logging.error("Registration of sample %s in database FAILED! Database is not updated.", sample_id)
        else:
            logging.info("User requested sample NOT to be registered in the database.")

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
