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

from lib.utils import read_fasta, get_the_data_from_the_input

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

    args.add_argument("--connstring",
                      "-c",
                      type=str,
                      metavar="CONNECTION",
                      required=True,
                      dest="db",
                      help="REQUIRED. Connection string for db.")

    args.add_argument("--reference",
                      type=str,
                      metavar="FASTAFILE",
                      required=True,
                      dest="reference",
                      help="""REQUIRED. Fasta reference file.""")

    args.add_argument("--refvar",
                      "-v",
                      metavar="JSONFILE",
                      required=True,
                      type=str,
                      dest="input",
                      help="""REQUIRED. Path to a Phenix-made json file containing variants
for the reference. Will contain Ns in positions that can't be mapped and these
will be globally ignored when calculating distances.""")

    args.add_argument("--ref-name",
                      "-r",
                      type=str,
                      metavar="NAME",
                      default=None,
                      dest="ref_name",
                      help="The name of the reference to go into the db [default: reference file name before 1st dot]")

    args.add_argument("--exclude",
                      "-e",
                      type=str,
                      metavar="BEDFILE",
                      default=None,
                      dest="exclude",
                      help="""Provide a bed file with regions to exclude globally. It needs to be a file with the
three required fields according to the specification (https://genome.ucsc.edu/FAQ/FAQformat.html#format1)
[default: don't exclude further regions]""")

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

        sql = "SELECT pk_id FROM samples"
        cur.execute(sql)
        if cur.rowcount > 0:
            logging.error("This is not an empty database.")
            return 1

        # open and read fasta reference file
        try:
            with open(args['reference'], 'r') as fa:
                dRefseq = read_fasta(fa)
        except IOError:
            logging.error("File not found: %s", args['reference'])
            return 1

        logging.info("%i contigs found in fasta reference.", len(dRefseq.keys()))

        # put contigs of reference into database
        contigs = {}
        for con in dRefseq.keys():
            sql = "INSERT INTO contigs (name, length) VALUES (%s, %s) RETURNING pk_id"
            cur.execute(sql, (con, len(dRefseq[con])))
            con_pkid = cur.fetchone()[0]
            contigs[con] = con_pkid

        # read the data from the json, fasta option does obviously not exist for ref
        args['format'] = 'json'
        data = get_the_data_from_the_input(args)
        if data == None:
            logging.error("An error occured getting the data from the input.")
            return 1

        # set refname if required
        if args['ref_name'] == None:
            args['ref_name'] = os.path.basename(args['reference']).split('.')[0]

        # get additional ignore postions from bed file if present
        exclude_regions = None
        if args['exclude'] != None:
            exclude_regions = get_exclude_regions_from_bed(args['exclude'], contigs)
            if exclude_regions == None:
                logging.error("An error occured with the exclude bed file.")
                return 1

        # ... make an entry in the samples table and get the primary sample id
        sql = "INSERT INTO samples (sample_name, date_added) VALUES (%s, %s) RETURNING pk_id"
        cur.execute(sql, (args['ref_name'], datetime.now(), ))
        ref_pkid = cur.fetchone()[0]

        logging.info("Created new sampe with id %s. ", ref_pkid)

        for con, condata in data['positions'].iteritems():
            # get the pk of this contig
            try:
                contig_pkid = contigs[con]
            except KeyError:
                logging.error("Contig %s which is in the json file was not found in the database. Does this sample belong in this database?", con)
                return 1

            # the reference can't have gaps, gaps in the vcf indicate a region where no reads have mapped back to
            # we treat these as Ns
            ref_ign_pos = set(condata['N']).union(set(condata['-']))
            # add additional global ignore positions as ref Ns to the database
            if exclude_regions != None and exclude_regions.has_key(con):
                ref_ign_pos.update(exclude_regions[con])

            # make one entry per contig in the variants table
            sql = "INSERT INTO variants (fk_sample_id, fk_contig_id, n_pos, a_pos, c_pos, g_pos, t_pos, gap_pos) VALUES (%s, %s, %s, %s, %s, %s, %s, %s)"
            cur.execute(sql, (ref_pkid, contig_pkid, list(ref_ign_pos), [], [], [], [], [], ))

            logging.info("Inserted for contig %s : Ns: %s",
                         contig_pkid,
                         len(ref_ign_pos))

        # per definition the reference is in cluster 1,1,1,1,1,1,1
        sql = "INSERT INTO sample_clusters (fk_sample_id, t0, t5, t10, t25, t50, t100, t250) VALUES (%s, 1, 1, 1, 1, 1, 1, 1)"
        cur.execute(sql, (ref_pkid, ))

        # and these are the stats for these clusters
        for t_lvl in ["t0", "t5", "t10", "t25", "t50", "t100", "t250"]:
            sql = "INSERT INTO cluster_stats (cluster_level, cluster_name, nof_members, nof_pairwise_dists) VALUES (%s, 1, 1, 0)"
            cur.execute(sql, (t_lvl, ))

        conn.commit()

    except psycopg2.Error as e:
         logging.error("Database reported error: %s" % (str(e)))
    finally:
        # close all dbs
        cur.close()
        conn.close()

    return 0

# --------------------------------------------------------------------------------------------------

def get_exclude_regions_from_bed(bedfile, contigs):
    """
    Parse a bed file and get the ignore regions as sets of each position.

    Parameters
    ----------
    bedfile: str
        the name of the bed file
    contigs: dict
        {contig_name: id, ...}

    Returns
    -------
    exclude_regions: dict
        {"contig1": set(1, 2, 3, 4, 5, 8, 9, 10, ...),
         "contig2": set(11, 12, 13, 14, 15, 18, 19, 110, 111, ...),
         ...}
    or None if problem
    """

    exclude_regions = {}
    try:
        with open(bedfile, 'r') as exbedf:
            for line in exbedf:
                cols = [x.strip() for x in line.strip().split('\t')]
                if len(cols) < 3:
                    logging.error("An error occured with the format of your exclude bed file.")
                    return None
                ex_contig = cols[0]
                try:
                    _ = contigs[ex_contig]
                except KeyError:
                    logging.error("A the contig name %s in the bed file does not correspond to the reference provided.",
                                  ex_contig)
                    return None
                try:
                    ex_start = int(cols[1]) + 1
                    ex_stop = int(cols[2])
                except ValueError:
                    logging.error("An error occured with the format of your exclude bed file.")
                    return None
                try:
                    exclude_regions[ex_contig].update(range(ex_start, ex_stop + 1))
                except KeyError:
                    exclude_regions[ex_contig] = set(range(ex_start, ex_stop + 1))

    except IOError:
        logging.error("An error occured reading from this file: %s", bedfile)
        return None

    return exclude_regions

# --------------------------------------------------------------------------------------------------

if __name__ == "__main__":
    exit(main(vars(get_args().parse_args())))
