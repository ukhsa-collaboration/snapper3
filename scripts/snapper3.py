#!/usr/bin/env python
'''
Single script for running other scripts in the project.

'''

import argparse
from argparse import RawTextHelpFormatter
import logging
import os

import add_sample
import add_reference
import cluster_sample
import get_alignment
import remove_sample
import get_closest
import export_sample_variants
import get_history
import make_tree
import update_db_trees

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

    subparsers.add_parser("add_reference",
                          description=add_reference.get_desc(),
                          help="Add the reference to an empty database.",
                          parents=[add_reference.get_args()],
                          formatter_class=RawTextHelpFormatter,
                          add_help=False)

    subparsers.add_parser("cluster_sample",
                          description=cluster_sample.get_desc(),
                          help="Put a sample into clusters.",
                          parents=[cluster_sample.get_args()],
                          formatter_class=RawTextHelpFormatter,
                          add_help=False)

    subparsers.add_parser("get_alignment",
                          description=get_alignment.get_desc(),
                          help="Get aligned sequences for a set of samples.",
                          parents=[get_alignment.get_args()],
                          formatter_class=RawTextHelpFormatter,
                          add_help=False)

    subparsers.add_parser("remove_sample",
                          description=remove_sample.get_desc(),
                          help="Remove a sample from the database.",
                          parents=[remove_sample.get_args()],
                          formatter_class=RawTextHelpFormatter,
                          add_help=False)

    subparsers.add_parser("get_closest",
                          description=get_closest.get_desc(),
                          help="Get the closest n samples from the database.",
                          parents=[get_closest.get_args()],
                          formatter_class=RawTextHelpFormatter,
                          add_help=False)

    subparsers.add_parser("export_sample_variants",
                          description=export_sample_variants.get_desc(),
                          help="Export variants for a given sample in json format.",
                          parents=[export_sample_variants.get_args()],
                          formatter_class=RawTextHelpFormatter,
                          add_help=False)

    subparsers.add_parser("get_history",
                          description=get_history.get_desc(),
                          help="Get a sample's SNP address history.",
                          parents=[get_history.get_args()],
                          formatter_class=RawTextHelpFormatter,
                          add_help=False)

    subparsers.add_parser("make_tree",
                          description=make_tree.get_desc(),
                          help="Make a phylogenetic tree for a group of samples.",
                          parents=[make_tree.get_args()],
                          formatter_class=RawTextHelpFormatter,
                          add_help=False)

    subparsers.add_parser("update_db_trees",
                          description=update_db_trees.get_desc(),
                          help="Update the trees stored in a database.",
                          parents=[update_db_trees.get_args()],
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
    elif args["cmd"] == "cluster_sample":
        return cluster_sample.main(args)
    elif args["cmd"] == "get_alignment":
        return get_alignment.main(args)
    elif args["cmd"] == "remove_sample":
        return remove_sample.main(args)
    elif args["cmd"] == "get_closest":
        return get_closest.main(args)
    elif args["cmd"] == "export_sample_variants":
        return export_sample_variants.main(args)
    elif args["cmd"] == "add_reference":
        return add_reference.main(args)
    elif args["cmd"] == "get_history":
        return get_history.main(args)
    elif args["cmd"] == "make_tree":
        return make_tree.main(args)
    elif args["cmd"] == "update_db_trees":
        return update_db_trees.main(args)
    return 1

# -------------------------------------------------------------------------------------------------

if __name__ == '__main__':
    exit(main())
