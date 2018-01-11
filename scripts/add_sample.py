import json
import logging
import gzip
import re
import os
from datetime import datetime

import argparse
from argparse import RawTextHelpFormatter
import psycopg2
from psycopg2.extras import DictCursor

from lib.utils import get_the_data_from_the_input, calculate_nless_n50

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
                      metavar="FILE",
                      required=True,
                      type=str,
                      dest="input",
                      help="REQUIRED. Path to a input file.")

    args.add_argument("--format",
                      "-f",
                      metavar="FORMAT",
                      required=True,
                      type=str,
                      dest="format",
                      choices=['json', 'fasta'],
                      help="REQUIRED. Choose from 'json' or 'fasta'.")

    args.add_argument("--connstring",
                      "-c",
                      type=str,
                      metavar="CONNECTION",
                      required=True,
                      dest="db",
                      help="REQUIRED. Connection string for db.")

    args.add_argument("--refname",
                      "-r",
                      type=str,
                      metavar="REFNAME",
                      required=True,
                      dest="refname",
                      help="REQUIRED. The sample_name of the reference genome in the database.")

    args.add_argument("--sample-name",
                      "-s",
                      type=str,
                      metavar="NAME",
                      default=None,
                      dest="sample_name",
                      help="The name of the sample to go into the db [default: input file name before 1st dot]")

    args.add_argument("--reference",
                      type=str,
                      metavar="FASTAFILE",
                      default=None,
                      dest="reference",
                      help="""Path to reference for this sample. Must be the same as used for the database.
REQUIRED when format is fasta, else ignored.""")

    args.add_argument("--min-coverage",
                      type=float,
                      metavar="FLOAT",
                      default=None,
                      dest="mcov",
                      help="""Minimum coverage required to allow sample in database. Only applicable with json
format, ignored for fasta.  This will check for the coverageMetaData annotation calculated by Phenix
and only allow the sample in, if the mean coverage is >= this value.
[default: do not check this]""")

    args.add_argument("--min-nless-n50",
                      type=int,
                      metavar="INT",
                      default=None,
                      dest="nless",
                      help="""Minimum N50 of N-less sequence required to allow sample in database. This will check
for the nlessnessMetaData annotation calculated by Phenix and only allow the sample in, if the n50 is >= this value.
If adding the sample from fasta it will be calculated on the fly. If adding from json and annotation is absent => error.
[default: do not check this]""")

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

    if args['format'] == 'fasta' and args['reference'] == None:
        logging.error("--reference option is REQUIRED when using fasta format input.")
        return 1

    data = get_the_data_from_the_input(args)

    if data == None:
        logging.error("An error occured getting the data from the input.")
        return 1

    # if format is json, check if there are any annotations to check and check them
    if args['format'] == 'json':
        if args['mcov'] != None:
            if data['annotations'].has_key('coverageMetaData') == False:
                logging.error("Was asked to check coverage but no coverage annotation found in json file.")
                return 1
            cov_info = dict(item.split("=") for item in data['annotations']['coverageMetaData'].split(","))
            logging.info("The mean coverage for this sample is: %s", cov_info['mean'])
            if float(cov_info['mean']) < args['mcov']:
                logging.error("The mean coverage for this sample is below the user specified threshold.")
                return 1
        if args['nless'] != None:
            if data['annotations'].has_key('nlessnessMetaData') == False:
                logging.error("Was asked to check nlessness but no nlessness annotation found in json file.")
                return 1
            nless_info = dict(item.split("=") for item in data['annotations']['nlessnessMetaData'].split(","))
            logging.info("The N-less N50 for this sample is: %s", nless_info['n50'])
            if int(nless_info['n50']) < args['nless']:
                logging.error("The N-less N50 for this sample is below the user specified threshold.")
                return 1

    # calculate nless annotation for fasta input if required
    if args['format'] == 'fasta' and args['nless'] != None:
        n50 = calculate_nless_n50(data, args['reference'])
        logging.info("The N-less N50 for this sample is: %s", n50)
        if n50 < args['nless']:
            logging.error("The mean coverage for this sample is below the user specified threshold.")
            return 1

    if args['sample_name'] == None:
        args['sample_name'] = os.path.basename(args['input']).split('.')[0]

    ngs_id = None
    molis_id = None
    # if this is the format <int>_H<int>-[12] add ngs_id and molis_id into db
    pat="^[0-9]+_H[0-9]+-[12]$"
    if re.search(pat, args['sample_name']) != None:
        ngs_id = int(args['sample_name'].split('_')[0])
        molis_id = args['sample_name'].split('_')[1][:-2]

    try:
        # open db
        conn = psycopg2.connect(args['db'])
        cur = conn.cursor(cursor_factory=psycopg2.extras.DictCursor)

        sql = "SELECT pk_id FROM samples WHERE sample_name=%s"
        cur.execute(sql, (args['sample_name'], ))
        if cur.rowcount > 0:
            logging.error("A sample with name %s already exists in the database.", args['sample_name'])
            return 1

        # get dict of contigs
        sql = "SELECT pk_id, name FROM contigs"
        cur.execute(sql)
        rows = cur.fetchall()
        contigs = {r['name']: r['pk_id'] for r in rows}

        # ... make an entry in the samples table and get the primary sample id
        sql = "INSERT INTO samples (sample_name, molis_id, ngs_id, date_added) VALUES (%s, %s, %s, %s) RETURNING pk_id"
        cur.execute(sql, (args['sample_name'], molis_id, ngs_id, datetime.now()))
        sample_pkid = cur.fetchone()[0]

        logging.info("Created new sampe with id %s. ", sample_pkid)

        for con, condata in data['positions'].iteritems():
            # get the pk of this contig
            try:
                contig_pkid = contigs[con]
            except KeyError:
                logging.error("Contig %s which is in the json file was not found in the database. Does this sample belong in this database?", con)
                return 1

            # get the positions on this contig ignored in the reference (n_pos) from the db
            sql = "SELECT v.n_pos FROM variants v, samples s WHERE v.fk_sample_id=s.pk_id AND s.sample_name=%s AND v.fk_contig_id=%s"
            cur.execute(sql, (args['refname'], contig_pkid, ))
            if cur.rowcount != 1:
                logging.error("Not exactly one row found in variants for sample %s and contig id %s.", args['refname'], contig_pkid)
                return 1
            res = cur.fetchone()[0]
            if res == None:
                ref_ign_pos = set()
            else:
                ref_ign_pos = set(res)

            # remove positions ignored in the reference from the variant positions registered for this sample
            # this means two things:
            # a) What is in the database is not an accurate representation of the sample's variants.
            # b) We never consider these regions when we calculate distances. <- This is what we want.
            a_pos = set(condata['A']) - ref_ign_pos
            c_pos = set(condata['C']) - ref_ign_pos
            g_pos = set(condata['G']) - ref_ign_pos
            t_pos = set(condata['T']) - ref_ign_pos
            n_pos = set(condata['N']) - ref_ign_pos
            gap_pos = set(condata['-']) - ref_ign_pos

            sql = "INSERT INTO variants (fk_sample_id, fk_contig_id, a_pos, c_pos, g_pos, t_pos, n_pos, gap_pos) VALUES (%s, %s, %s, %s, %s, %s, %s, %s)"
            cur.execute(sql, (sample_pkid, contig_pkid, list(a_pos), list(c_pos), list(g_pos), list(t_pos), list(n_pos), list(gap_pos), ))

            logging.info("Inserted for contig %s : As: %s, Cs: %s:, Gs: %s, Ts: %s, Ns: %s, gaps: %s",
                         contig_pkid,
                         len(a_pos),
                         len(c_pos),
                         len(g_pos),
                         len(t_pos),
                         len(n_pos),
                         len(gap_pos))

        # missing contigs are contigs that are all ref, which might be absent from the json file
        # if there are any like that, we need to add empty rows for them
        missing_contigs = set(contigs.keys()).difference(data['positions'].keys())
        if len(missing_contigs) > 0:
            for mc in missing_contigs:
                sql = "INSERT INTO variants (fk_sample_id, fk_contig_id, a_pos, c_pos, g_pos, t_pos, n_pos, gap_pos) VALUES (%s, %s, %s, %s, %s, %s, %s, %s)"
                cur.execute(sql, (sample_pkid, contigs[mc], [], [], [], [], [], [], ))
                logging.info("Inserted for contig %s : As: %i, Cs: %i:, Gs: %i, Ts: %i, Ns: %i, gaps: %i", contigs[mc], 0, 0, 0, 0, 0, 0)

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
