"""
File contains some stuff to make alignments in snapper3.

author: ulf.schaefer@phe.gov.uk

"""

import logging
import os
import sys

from lib.utils import read_fasta

import psycopg2
from psycopg2.extras import DictCursor

# --------------------------------------------------------------------------------------------------

def add_reference_data(ref, all_contig_data):
    """
    Since we don't have it in the database, add the refbases for all relevant positions from
    the reference fasta file.

    Parameters
    ----------
    ref: dict
        fasta file in dict as read by readfasta
    all_contig_data: dict
        see get_alignment.main for what it looks like
    Returns
    -------
    0 if fine
    None if fail
    """

    if sorted(ref.keys()) != sorted(all_contig_data.keys()):
        logging.error("Contig names don't match between the database and the file passed in the --reference parameter.")
        return None

    for contig, data in all_contig_data.iteritems():
        data['reference'] = {'A': set(), 'C': set(), 'G': set(), 'T': set(), 'N': set(), '-': set()}
        for sam in data.keys():
            for n in data[sam].keys():
                for x in data[sam][n]:
                    refbase = ref[contig][x-1].upper()
                    data['reference'][refbase].add(x)

    return 0

# --------------------------------------------------------------------------------------------------

def get_data_from_db(db, samples_in, refname):
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

        sql  = "SELECT pk_id FROM samples WHERE sample_name=%s"
        cur.execute(sql, (refname, ))
        if cur.rowcount != 1:
            logging.error("Reference name not found in database.")
            return None
        refid = cur.fetchone()[0]

        print refid

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

            # get the positions on this contig ignored in the reference (n_pos) from the db
            sql = "SELECT n_pos FROM variants WHERE fk_sample_id=%s AND fk_contig_id=%s"
            cur.execute(sql, (refid, con_id, ))
            if cur.rowcount != 1:
                logging.error("Not exactly one row found in variants for reference on contig id %s.", con_id)
                return None
            res = cur.fetchone()[0]
            if res == None:
                ref_ign_pos = set()
            else:
                ref_ign_pos = set(res)

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
                all_contig_data[con_name][sam_name]['-'] = set(r['gap_pos'])
                all_contig_data[con_name][sam_name]['N'] = set(r['n_pos'])
                all_contig_data[con_name][sam_name]['N'].update(ref_ign_pos) # add reference ignore positions back in

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

def process_bed_file(args, all_contig_data):
    """
    Process a bed file with position intervals to include or exclude.

    Parameters:
    -----------
    args: dict
        arguments dictionary
    all_contig_data: dict
        all data as described abobe in main

    Returns:
    --------
    0, but all_contig_data is updated
    """

    all_samples = set()
    for contig in all_contig_data.keys():
        all_samples.update(all_contig_data[contig].keys())

    incl = True if args['include'] else False

    with open(args["exclude"]) as fp:
        for line in fp:
            data = line.strip().split("\t")
            try:
                for sam in all_samples:
                    for nuc in all_contig_data[data[0]][sam].keys():
                        if incl == True:
                            all_contig_data[data[0]][sam][nuc].intersection_update(range(int(data[1]),
                                                                                         int(data[2]) + 1))
                        else:
                            all_contig_data[data[0]][sam][nuc].difference_update(range(int(data[1]),
                                                                                       int(data[2]) + 1))
            except KeyError:
                logging.error("Wrong contig in bed file: %s. Ignoring.", data[0])

    return 0

# --------------------------------------------------------------------------------------------------

def remove_columns(t, character, all_contig_data):
    """
    Remove columns form the alignment if the character [N, gap] is above
    the fraction t in a column.

    Parameters:
    -----------
    t: float
        threshold between 0.0 and 1.0
    character: str
        'N' or '-'
    all_contig_data: dict
        all data as described abobe in main

    Returns:
    --------
    0, but all_contig_data is updated
    """

    for (contig, data) in all_contig_data.iteritems():
        all_pos = set()
        for nuc in data['reference']:
            all_pos.update(data['reference'][nuc])
        # number of samples not considering the reference
        nof_samples = len(data.keys()) - 1
        to_remove = set()
        for pos in all_pos:
            ns = 0
            for sam in data.keys():
                if sam != 'reference':
                    try:
                        if pos in data[sam][character]:
                            ns += 1
                    except KeyError:
                        # happens when sam doesn't have any Ns
                        pass
            if float(ns) / nof_samples > t:
                to_remove.add(pos)

        # remove postions
        if len(to_remove) > 0:
            for sam in data.keys():
                for nuc in data[sam].keys():
                    data[sam][nuc].difference_update(to_remove)

    return 0

# --------------------------------------------------------------------------------------------------

def remove_samples(args, option, character, all_contig_data):
    """
    Remove samples form the alignment if the character [N, gap] is above
    a given fraction in in samples (relative to the size of the genome.)

    Parameters:
    -----------
    args: dict
        arguments dictionary
    option: str
        'sample_Ns' or 'sample_gaps'
    character: str
        'N' or '-'
    all_contig_data: dict
        all data as described abobe in main

    Returns:
    --------
    0, but all_contig_data is updated
    """

    # get the alignment length
    align_len = args["reflength"]

    # sum up the number of Ns or gaps for each sample
    ns_per_sample = {}
    for (contig, data) in all_contig_data.iteritems():
        for sam in data.keys():
            if ns_per_sample.has_key(sam) == False:
                ns_per_sample[sam] = 0
            try:
                ns_per_sample[sam] += len(data[sam][character])
            except KeyError:
                # this happens when sam has no Ns
                pass

    # calculate proportion of Ns or gaps
    for sam in ns_per_sample.keys():
        ns_per_sample[sam] = float(ns_per_sample[sam]) / align_len

    # if this is set to auto calculate the threshold from the mean + 2 times the stddev
    t  = 0.0
    if args[option] == 'auto':
        m = sum(ns_per_sample.values()) / len(ns_per_sample)
        ssd = sum((x-m)**2 for x in ns_per_sample.values())
        variance = ssd / len(ns_per_sample)
        sd = sqrt(variance)
        t = m + (args['sample_Ns_gaps_auto_factor']*sd)
    else:
        t = args[option]

    # remove all samples thath are bigger than the thresdold
    removals = False
    for sam in ns_per_sample.keys():
        if ns_per_sample[sam] > t:
            for (contig, data) in all_contig_data.iteritems():
                logging.info("Removing sample %s, because it has %.3f %ss", sam, ns_per_sample[sam], character)
                del data[sam]
                removals = True

    # tidy up
    if removals == True:
        for (contig, data) in all_contig_data.iteritems():
            # get all positions in the reference and all positions in all other samples
            ref_pos = set()
            var_pos = set()
            for sam in data.keys():
                if sam == 'reference':
                    for nuc in data[sam].keys():
                        ref_pos.update(data[sam][nuc])
                else:
                    for nuc in data[sam].keys():
                        var_pos.update(data[sam][nuc])
            # get all positions that are only in the reference
            only_ref_pos = ref_pos.difference(var_pos)
            # remove the positions only in the reference from the reference
            # because we probably removed the sample with a variant at those positions
            for nuc in data['reference'].keys():
                data['reference'][nuc].difference_update(only_ref_pos)

    return 0

# -------------------------------------------------------------------------------------------------

def output_per_sample_stats(all_contig_data):
    """
    Output stats per sample in the aligment.

    Parameters:
    -----------
    all_contig_data: dict
        all data as described abobe in main

    Returns:
    --------
    0
    """

    all_contigs = all_contig_data.keys()
    all_samples = set()
    for contig in all_contigs:
        all_samples.update(all_contig_data[contig].keys())
    for sam in all_samples:
        tot = 0
        ns = 0
        gaps = 0
        mut = 0
        mix = 0
        for contig in all_contigs:
            for k in all_contig_data[contig][sam].keys():
                tot += len(all_contig_data[contig][sam][k])
            try:
                ns += len(all_contig_data[contig][sam]['N'])
            except KeyError:
                pass
            try:
                gaps += len(all_contig_data[contig][sam]['-'])
            except KeyError:
                pass
            for k in ['A', 'C', 'G', 'T']:
                try:
                    mut += len(all_contig_data[contig][sam][k])
                except KeyError:
                    pass
            mixbases = set(all_contig_data[contig][sam].keys()).difference(['A', 'C', 'G', 'T', 'N', '-'])
            for m in mixbases:
                mix += len(all_contig_data[contig][sam][m])
        if sam == 'reference':
            mut = 0
        sys.stdout.write("%s\tN: %i, mut: %i, mix: %i, gap: %i, total: %i\n" %(sam, ns, mut, mix, gaps, tot))

    return 0

# --------------------------------------------------------------------------------------------------
