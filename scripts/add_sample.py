import json
import logging
import argparse
from argparse import RawTextHelpFormatter
import psycopg2
from psycopg2.extras import DictCursor

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

    data = None
    try:
        with open(args['input']) as data_file:
            try:
                data = json.load(data_file)
            except ValueError:
                logging.error("Data in %s is corrupted.", args['input'])
                return 1
    except IOError:
        logging.error("Could not open file %s", args['input'])
        return 1

    print data

    try:
        # open db
        conn = psycopg2.connect(args['db'])
        cur = conn.cursor(cursor_factory=psycopg2.extras.DictCursor)

        # put code here


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
