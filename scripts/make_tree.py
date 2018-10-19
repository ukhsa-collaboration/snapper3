import logging
import os
import argparse
from argparse import RawTextHelpFormatter

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

    return r'''For a given set of samples, this can either return a NJ or an ML tree.'''

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

    args.add_argument("--method",
                      "-m",
                      type=str,
                      metavar="TREEMETHOD",
                      required=True,
                      choices=["NJ", "ML"],
                      dest="method",
                      help="""The phylogenetic method for making the tree.
Choose from 'ML' or 'NJ'. REQUIRED.""")

    args.add_argument("--output",
                      "-o",
                      type=str,
                      metavar="FILENAME",
                      default='tree.nwk',
                      dest="output",
                      help="""The name of the output file. [Default: tree.nwk]""")

    args.add_argument("--distance-matrix",
                      "-d",
                      type=str,
                      metavar="FILENAME",
                      default=None,
                      dest="dm",
                      help="""Store distance matrix as csv in this file.
Ignored for ML. [Default: Do not store]""")

    args.add_argument("--samples",
                      "-s",
                      type=str,
                      metavar="SAMPLES",
                      default=None,
                      dest="samples",
                      help="""List of samples to put in the tree. Can be either file
or a comma-separated list w/o blanks. [Default: None, but then --clusters has to be used.]""")

    args.add_argument("--clusters",
                      "-l",
                      type=str,
                      metavar="CLUSTERS",
                      default=None,
                      dest="clusters",
                      help="""List of clusters to put in the tree. List of key:value pairs
, e.g. t100:25,t50:13,t50:14 [Default: None, but then --samples has to be used.]""")

    args.add_argument("--reference",
                      "-r",
                      type=str,
                      default=None,
                      dest="ref",
                      help="Path to reference specified (FASTA). Required for ML, else ignored.")

    args.add_argument("--refname",
                      type=str,
                      default=None,
                      help="""The name of the reference in the database.
Required for ML, else ignored.""")

    args.add_argument("--remove-ref",
                      action='store_true',
                      help="""Remove the reference from the ML tree [Default: reference is always
in ML tree]. Ignored for NJ. NJ have the reference in
when it's specified in --clusters or --samples.""")

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

    if args['samples'] == None and args['clusters'] == None:
        logging.error("Either --samples or --clusters need to be specified.")
        return 1

    if args['method'] == 'ML':
        if args['ref'] == None or args['refname'] == None:
            logging.error("For ML trees, both --reference and --refname are required.")
            return 1
        if args['dm'] != None:
            logging.warning("Distance matrix is not going to be saved for ML trees.")

    samples = None
    if args['samples'] != None:
        if os.path.exists(args['samples']) == True and os.path.isfile(args['samples']) == True:
            logging.info("Reading sample names from file: %s", args['samples'])
            with open(args['samples'], 'r') as f:
                samples = [line.strip() for line in f]
        else:
            logging.info("Trying to turns this in a list of sample names: %s", args['samples'])
            samples = args['samples'].split(',')

        logging.info("Using the following list of samples names: %s", samples)

    clusters = {}
    if args['clusters'] != None:
        for kv_pair in args['clusters'].split(","):
            pair = kv_pair.split(":")
            if len(pair) != 2:
                logging.error("Clusters should be separated by ':'. Please consider: %s", str(pair))
                return 1
            try:
                clusters[pair[0]].append(pair[1])
            except KeyError:
                clusters[pair[0]] = [pair[1]]

        logging.info("Using the following lists of clusters: %s", str(clusters))
    else:
        clusters = None

    result = None
    with SnapperDBInterrogation(conn_string=args['db'], fusion=args['fusion']) as sdbi:
        try:
            result = sdbi.get_tree(samples,
                                   clusters,
                                   args['method'],
                                   ref=args['ref'],
                                   refname=args['refname'],
                                   dm=args['dm'],
                                   rmref=args['remove_ref'])
            with open(args['output'], 'w') as f:
                f.write(result)
                logging.info("Tree written to file: %s", args['output'])
        except SnapperDBInterrogationError as e:
            logging.error(e)
        else:
            logging.info("Completed successfully.")

    return 0

# end of main --------------------------------------------------------------------------------------

if __name__ == "__main__":
    exit(main(vars(get_args().parse_args())))
