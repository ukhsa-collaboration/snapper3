import logging

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

    return r'''This connects to the database and returns the closest N samples to a given input sample.'''

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
                      help="Name of sample to query. REQUIRED.")

    args.add_argument("--neighbours",
                      "-n",
                      type=int,
                      default=1,
                      dest="neighbours",
                      help="""The number of neighbours to return. Additional samples that have
the same distance as the nth sample are also returned. [DEFAULT: 1]""")

    args.add_argument("--distance",
                      "-d",
                      type=int,
                      default=None,
                      dest="distance",
                      help="""Return all samples closer or equal to this distance. [DEFAULT: None, i.e.
return the n closest samples specified by -n option.]""")

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

    with SnapperDBInterrogation(conn_string=args['db']) as sdbi:
        try:

            result = None
            if args['distance'] == None:
                result = sdbi.get_closest_samples(args['sample'], args['neighbours'])
            else:
                result = sdbi.get_samples_below_threshold(args['sample'], args['distance'])

            for i, x in enumerate(result):
                print i+1, x

        except SnapperDBInterrogationError as e:
            logging.error(e)

    return 0

# end of main --------------------------------------------------------------------------------------

if __name__ == "__main__":
    exit(main(vars(get_args().parse_args())))
