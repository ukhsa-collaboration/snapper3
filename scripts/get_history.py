import logging
import sys
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

    result = None
    with SnapperDBInterrogation(conn_string=args['db']) as sdbi:
        try:
            result = sdbi.get_sample_history(args['sample'])

            sys.stdout.write("Current SNP address: %s\n" % (result['current_snad']))
            if len(result['history']) == 0:
                sys.stdout.write("SNP address has never changed.\n")
            else:
                for hi in sorted(result['history'], key=lambda k: k['time'], reverse=True):
                    sys.stdout.write("SNP address change @ %s: %s -> %s\n" \
                                     % (hi['time'].strftime("%Y-%m-%d %H:%M:%S"), hi['old'], hi['new']))

        except SnapperDBInterrogationError as e:
            logging.error(e)

    return 0

# end of main --------------------------------------------------------------------------------------

if __name__ == "__main__":
    exit(main(vars(get_args().parse_args())))
