import logging

import requests
import argparse
from argparse import RawTextHelpFormatter
import psycopg2
from psycopg2.extras import DictCursor

import lib.snapperdb as sndb
import lib.registration as regis
import lib.merging as merging
from lib.distances import get_all_pw_dists, get_relevant_distances
from lib.distances import get_distances_precalc, get_distances_fusion

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

    args.add_argument("--precalc-distances",
                      "-p",
                      type=str,
                      metavar="JSONFILE",
                      default=None,
                      dest="precalc",
                      help="""Json file with precalculated distances.
[Default: None => calculate all distances on the db server now]""")

    args.add_argument("--no-zscore-check",
                      action='store_true',
                      help="""Do not perform checks and just add the sample. It's fine.
[Default: Perform checks.]""")

    args.add_argument("--with-registration",
                      action='store_true',
                      help="""Register the clustering for this sample in the database
and update the cluster stats. [Default: Do not register.]""")

    args.add_argument("--force-merge",
                      action='store_true',
                      help="""Add the sample even if it causes clusters to merge.
[Default: Do not add if merge required.]""")

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
            logging.error("Could not find a cluster-able samples with this name.")
            return 1

        rc, snadd = sndb.check_duplicate_clustering(cur, sample_id)
        if rc != 0:
            logging.error("Sample %s with sample id %s was already clustered with snp address: %s",
                          args['sample_name'],
                          sample_id,
                          snadd)
            return 1

        if args['precalc'] != None and args['fusion'] != None:
            logging.error("Parameters --precalc-distances and --fusion are mutually exclusive. Use only one of them (or neither).")
            return 1

        logging.info("Processing sample %s with id %i", args['sample_name'], sample_id)
        logging.info("Calculating distances to all other samples now. Patience!")

        # calculate the relevant distances
        if args['precalc'] != None:
            distances = get_distances_precalc(cur, sample_id, args['sample_name'], args['precalc'])
        elif args['fusion'] != None:
            distances = get_distances_fusion(cur, sample_id, args['fusion'])
        else:
            distances = get_relevant_distances(cur, sample_id)

        if distances == None:
            logging.error("Could not get distances. :-(")
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
                      new_snad[6],
                      new_snad[5],
                      new_snad[4],
                      new_snad[3],
                      new_snad[2],
                      new_snad[1],
                      new_snad[0])

        merges = merging.check_merging_needed(cur, distances, new_snad)

        logging.info("Merges that would be required to make this assignment: %s",
                     str([str(m) for m in merges.values()]))

        if args['force_merge'] == False and len(merges.keys()) > 0:
            logging.info("Exiting becaues there are merges. Database is not updated. Use --force-merge to add the sample (after checking it).")
            return 1

        if args['no_zscore_check'] == False:
            zscore_fail, zscore_info = sndb.check_zscores(cur, distances, new_snad, merges)

            if zscore_fail == None:
                logging.error("Could not calculate z-scores. :-(")
                return 1

            if zscore_fail == True:
                for s in zscore_info:
                    logging.info(s)
                logging.error("z-score check for this assignment has failed. Database is not updated.")
                return 1

            logging.info("All z-score checks passed for this assignment.")
        else:
            logging.info("User disabled z-score checks for this assignment.")

        logging.debug("Merges that would be required to make this assignment: %s",
                      str([str(m) for m in merges.values()]))

        if args['with_registration'] == True:

            levels = [0, 5, 10, 25, 50, 100, 250]
            for lvl in merges.keys():
                merging.do_the_merge(cur, merges[lvl])
                # If merging cluster a and b, the final name of the merged cluster can be either
                # a or b. So we need to make sure the cluster gets registered into the final name
                # of the cluster and not into the cluster that has been deleted in the merge
                # operation.
                new_snad[levels.index(lvl)] = merges[lvl].final_name

            final_snad = regis.register_sample(cur,
                                               sample_id,
                                               distances,
                                               new_snad,
                                               args['no_zscore_check'])

            if final_snad != None:
                logging.info("Sample %s with sample_id %s was registered in the database with SNP address: %s-%s-%s-%s-%s-%s-%s",
                             args['sample_name'],
                             sample_id,
                             final_snad[6],
                             final_snad[5],
                             final_snad[4],
                             final_snad[3],
                             final_snad[2],
                             final_snad[1],
                             final_snad[0])

                # if we did not do zscore checks, because the sample would fail, we need to ignore it in the future
                if args['no_zscore_check'] == True:
                    sql = "UPDATE samples SET ignore_zscore=TRUE WHERE pk_id=%s"
                    cur.execute(sql, (sample_id, ))

            else:
                logging.error("Registration of sample %s in database FAILED! Database is not updated.",
                              sample_id)
                return 1
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
