import json
import logging
import gzip
import re
from datetime import datetime

import argparse
from argparse import RawTextHelpFormatter
import psycopg2
from psycopg2.extras import DictCursor

from lib.utils import check_json_format

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

    return r'''Takes variants for a sample in json format and adds them to the database.'''

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

    args.add_argument("--input",
                      "-i",
                      metavar="JSONFILE",
                      required=True,
                      type=str,
                      dest="input",
                      help="path to a json file")

    args.add_argument("--connstring",
                      "-c",
                      type=str,
                      metavar="CONNECTION",
                      required=True,
                      dest="db",
                      help="Connection string for db")

    args.add_argument("--sample-name",
                      "-s",
                      type=str,
                      metavar="NAME",
                      default=None,
                      dest="sample_name",
                      help="The name of the sample to go into the db [default: input file name before 1st dot]")


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

    data = None
    open_func = gzip.open if args['input'].endswith('.gz') == True else open
    try:
        with open_func(args['input']) as data_file:
            try:
                data = json.load(data_file)
            except ValueError:
                logging.error("Data in %s is corrupted.", args['input'])
                return 1
    except IOError:
        logging.error("Could not open file %s", args['input'])
        return 1

    if check_json_format(data) == False:
        logging.error("Data in %s is not in the correct format. Pleas use the latest version of Phenix to make this file from a filtered vcf.", args['input'])
        return 1

    if args['sample_name'] == None:
        args['sample_name'] = args['input'].split('.')[0]

    ngs_id = None
    molis_id = None

    # if this is the format <int>_H<int>-[12] add ngs_id and molis_id into db
    pat="^[0-9]+_H[0-9]+-[12]"
    if re.search(pat, args['sample_name']) != None:
        ngs_id = int(args['sample_name'].split('_')[0])
        molis_id = args['sample_name'].split('_')[1][:-2]

    print args
    print ngs_id
    print molis_id



    try:
        # open db
        conn = psycopg2.connect(args['db'])
        cur = conn.cursor(cursor_factory=psycopg2.extras.DictCursor)

        # put code here

        # ... make an entry in the samples table and get the primary sample id
        sql = "INSERT INTO samples (sample_name, molis_id, ngs_id, date_added) VALUES (%s, %s, %s, %s) RETURNING pk_id"
        cur.execute(sql, (args['sample_name'], molis_id, ngs_id, datetime.now()))
        sample_pkid = cur.fetchone()[0]

        print sample_pkid

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
