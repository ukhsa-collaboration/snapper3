#!/usr/bin/env python
'''
Single script for running other scripts in the project.

'''

import argparse
from argparse import RawTextHelpFormatter
import logging
import os

import add_sample

# -------------------------------------------------------------------------------------------------

def get_version():
    """
    Get the version from the VERSION file that should be one level up.
    Return 'n/a' if the file cannot be found

    Parameters
    ----------
    no inputs

    Returns
    -------
    version: str
        the version number
    """

    version_file = os.path.join(os.path.abspath(os.path.dirname(__file__)), "..", "VERSION")
    version = "N/A"

    if os.path.exists(version_file):
        try:
            with open(version_file) as fp:
                version = fp.next().strip()
        except IOError:
            pass
    return version

# -------------------------------------------------------------------------------------------------

def get_args():
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

    args = argparse.ArgumentParser()

    args.add_argument("--debug",
                      action="store_true",
                      help="More verbose logging (default: turned off).")

    args.add_argument("--version",
                      action="version",
                      version=get_version())

    subparsers = args.add_subparsers(dest='cmd')

    subparsers.add_parser("add_sample",
                          description=add_sample.get_desc(),
                          help="Add a sample to the database.",
                          parents=[add_sample.get_args()],
                          formatter_class=RawTextHelpFormatter,
                          add_help=False)

    return args

# -------------------------------------------------------------------------------------------------

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

    version = get_version()

    args = vars(get_args().parse_args())

    log_level = logging.DEBUG if args["debug"] else logging.INFO
    logging.basicConfig(format="[%(asctime)s] %(levelname)s: %(message)s",
                            level=log_level)

    logging.info("Version: %s", version)

    args["version"] = version

    if args["cmd"] == "add_sample":
        return add_sample.main(args)
    return 1

# -------------------------------------------------------------------------------------------------

if __name__ == '__main__':
    exit(main())
