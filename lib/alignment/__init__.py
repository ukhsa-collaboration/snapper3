"""
File contains some stuff to make alignments in snapper3.

author: ulf.schaefer@phe.gov.uk

"""

import logging
import os

from lib.utils import read_fasta

import psycopg2
from psycopg2.extras import DictCursor

# --------------------------------------------------------------------------------------------------

def add_reference_data(refname, all_contig_data):
    """
    Get the sets of variant positions from the database for all samples and all contigs.

    Parameters
    ----------
    db: str
        database conncetion string
    samples_in: list
        input list of sample names

    Returns
    -------
    all_contig_data: dict
        see get_alignment.main for what it looks like
    """

    if os.path.exists(refname) == False:
        logging.error("File not found: %s", refname)
        return None

    ref = {}
    with open(refname, 'r') as ref_file:
        ref = read_fasta(ref_file)

    if sorted(ref.keys()) != sorted(all_contig_data.keys()):
        logging.error("Contig names don't match between the database and the file passed in the --reference parameter.")
        return None

    for contig, data in all_contig_data.iteritems():
        data['reference'] = {'A': set(), 'C': set(), 'G': set(), 'T': set(), 'N': set(), '-': set()}
        for sam in data.keys():
            for n in data[sam].keys():
                for x in data[sam][n]:
                    refbase = ref[contig][x-1].upper()
                    data[reference][refbase].add(x)

    return 0

# --------------------------------------------------------------------------------------------------

def get_data_from_db(db, samples_in):
    """
    Get the sets of variant positions from the database for all samples and all contigs.

    Parameters
    ----------
    db: str
        database conncetion string
    samples_in: list
        input list of sample names

    Returns
    -------
    all_contig_data: dict
        see get_alignment.main for what it looks like
    """

    all_contig_data = {}

    try:
        # open db
        conn = psycopg2.connect(db)
        cur = conn.cursor(cursor_factory=psycopg2.extras.DictCursor)

        sql = "SELECT pk_id, sample_name FROM samples WHERE sample_name IN %s"
        cur.execute(sql, (tuple(samples_in), ))
        rows = cur.fetchall()
        samples = {r['pk_id']: r['sample_name'] for r in rows}

        miss = set(samples_in).difference(set(samples.values()))

        if len(miss) == len(samples_in):
            logging.error("None of your input samples could be found in the database.")
            return None
        if len(miss) > 0:
            logging.warning("These samples could not be found in the database: %s", str(list(miss)))

        sql = "SELECT pk_id, name FROM contigs"
        cur.execute(sql)
        rows = cur.fetchall()
        contigs = {r['pk_id']: r['name'] for r in rows}

        for con_id, con_name in contigs.iteritems():

            all_contig_data[con_name] = {}

            sql = "SELECT fk_sample_id, a_pos, c_pos, g_pos, t_pos, n_pos, gap_pos FROM variants WHERE fk_sample_id IN %s AND fk_contig_id=%s"
            cur.execute(sql, (tuple(samples.keys()), con_id, ))
            rows = cur.fetchall()
            for r in rows:
                sam_name = samples[r['fk_sample_id']]
                all_contig_data[con_name][sam_name] = {}
                all_contig_data[con_name][sam_name]['A'] = set(r['a_pos'])
                all_contig_data[con_name][sam_name]['C'] = set(r['c_pos'])
                all_contig_data[con_name][sam_name]['G'] = set(r['g_pos'])
                all_contig_data[con_name][sam_name]['T'] = set(r['t_pos'])
                all_contig_data[con_name][sam_name]['N'] = set(r['n_pos'])
                all_contig_data[con_name][sam_name]['-'] = set(r['gap_pos'])

        conn.commit()

    except psycopg2.Error as e:
         logging.error("Database reported error: %s" % (str(e)))
         return None
    finally:
        # close all dbs
        cur.close()
        conn.close()

    return all_contig_data

# --------------------------------------------------------------------------------------------------
