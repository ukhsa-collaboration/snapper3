import json
import logging
import gzip
from datetime import datetime

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

    return r'''Exports the variants for a sample in json format.'''

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

    args.add_argument("--sample-name",
                      "-s",
                      type=str,
                      metavar="NAME",
                      default=None,
                      dest="sample_name",
                      help="The name of the sample in the db. REQUIRED")

    args.add_argument("--connstring",
                      "-c",
                      type=str,
                      metavar="CONNECTION",
                      required=True,
                      dest="db",
                      help="Connection string for db. REQUIRED")

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

    try:
        # open db
        conn = psycopg2.connect(args['db'])
        cur = conn.cursor(cursor_factory=psycopg2.extras.DictCursor)

        # check sampleis in the database
        sql = "SELECT pk_id FROM samples WHERE sample_name=%s"
        cur.execute(sql, (args['sample_name'], ))
        if cur.rowcount != 1:
            logging.error("Could not get sample %s from database.", args['sample_name'])
            return 1
        sam_id = cur.fetchone()[0]

        # get the db name from the connstring for the annotation
        bla=args['db'][args['db'].find("dbname")+8:]
        dbname = bla[:bla.find("'")]

        data = {}
        data['annotations'] = {"dbname": dbname}
        data['positions'] = {}

        # get all contigs from the database
        sql = "SELECT pk_id, name FROM contigs"
        cur.execute(sql)
        rows = cur.fetchall()
        # for all contigs get variants info and add to data dict
        for r in rows:
            sql = "SELECT a_pos, c_pos, g_pos, t_pos, n_pos, gap_pos FROM variants WHERE fk_sample_id=%s AND fk_contig_id=%s"
            cur.execute(sql, (sam_id, r['pk_id'], ))
            if cur.rowcount != 1:
                logging.error("Could not get variants from database for sample %s on contig %s.", args['sample_name'], r['pk_id'])
                return 1

            var = cur.fetchone()
            data['positions'][r['name']] = {"G": var['g_pos'],
                                            "A": var['a_pos'],
                                            "T": var['t_pos'],
                                            "C": var['c_pos'],
                                            "N": var['n_pos'],
                                            "-": var['gap_pos']}

        # create json string and write to file
        filename = "%s.json.gz" % (args['sample_name'])
        logging.info("Writing data to file: %s", filename)
        json_string = json.dumps(data)
        with gzip.open(filename, mode='wb') as gzip_obj:
            gzip_obj.write(json_string)

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
