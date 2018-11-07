#!/usr/bin/env python

import sys
import os
import argparse
import logging
import psycopg2
from psycopg2.extras import DictCursor
from math import sqrt

from time import time

__version__ = '0.1'
__date__ = '30Sep2016'
__author__ = 'ulf.schaefer@phe.gov.uk'

# --------------------------------------------------------------------------------------------------

def parse_args():
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

    sDescription = 'version %s, date %s, author %s' %(__version__, __date__, __author__)

    parser = argparse.ArgumentParser(description=sDescription)

    parser.add_argument("--reference",
                         "-r",
                         type=str,
                         metavar="REFNAME",
                         required=True,
                         dest="ref",
                         help="The sample_name of the reference genome in the database.")

    parser.add_argument("--oldconnstring",
                        "-o",
                        type=str,
                        metavar="CONNECTION",
                        required=True,
                        dest="olddb",
                        help="Connection string for old db ('source')")

    parser.add_argument("--newconnstring",
                        "-n",
                        type=str,
                        metavar="CONNECTION",
                        required=True,
                        dest="newdb",
                        help="Connection string for new db ('target')")

    oArgs = parser.parse_args()
    return oArgs

# --------------------------------------------------------------------------------------------------

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
    oArgs = parse_args()

    logging.basicConfig(format="[%(asctime)s] %(levelname)s: %(message)s", level=logging.INFO)

    try:
        # open source db
        source_conn = psycopg2.connect(oArgs.olddb)
        source_cur = source_conn.cursor(cursor_factory=psycopg2.extras.DictCursor)

        # open target db
        target_conn = psycopg2.connect(oArgs.newdb)
        target_cur = target_conn.cursor(cursor_factory=psycopg2.extras.DictCursor)

        logging.info("Starting to migrate samples and clusters.")
        all_sample_names = migrate_samples_and_clusters(source_cur, target_cur)
        logging.info("Completed migration of samples and clusters.")

        if oArgs.ref in all_sample_names == False:
            logging.error("Reference name provided (%s) not found in the database.", oArgs.ref)
            raise SystemExit("critical error")

        all_contigs = migrate_contigs(source_cur, target_cur)

        logging.info("Starting to migrate variants.")
        migrate_variants(oArgs.ref, source_cur, target_cur, all_sample_names)
        logging.info("Completed variants migration.")

        logging.info("Reading distance matrix into memory.")
        dm = {}
        read_distance_matrix(source_cur, dm)
        logging.info("Completed reading distance matrix.")

        logging.info("Calculating per cluster stats.")
        calculate_per_cluster_stats(target_cur, dm)
        logging.info("Completed per cluster stats")

        logging.info("Calculating per sample cluster stats.")
        calculate_per_sample_cluster_stats(target_cur, dm)
        logging.info("Completed per sample cluster stats")

        logging.info("MIGRATION SUCCESSFULLY COMPLETED.")

        target_conn.commit()

    except SystemExit as e:
        logging.error("Could not complete migration because: %s" % (str(e)))
    except psycopg2.Error as e:
         logging.error("Database reported error: %s" % (str(e)))
    finally:
        # close all dbs
        source_cur.close()
        source_conn.close()
        target_cur.close()
        target_conn.close()

    return 0

# end of main --------------------------------------------------------------------------------------

def calculate_per_sample_cluster_stats(target_cur, dm):
    """
    For each sample in each of its clusters calculate the mean distance of this sample to all other
    cluster members and update the sample_clusters table.

    Parameters
    ----------
    target_cur: obj
        database cursor
    dm: dict
        container for distance matrix

    Returns
    -------
    always 0
    """

    # get all samples and their clusters from db
    sql = "SELECT c.pk_id, c.fk_sample_id, c.t0, c.t5, c.t10, c.t25, c.t50, c.t100, c.t250 FROM sample_clusters c, samples s WHERE s.ignore_sample=false AND s.ignore_zscore=false AND c.fk_sample_id=s.pk_id"
    target_cur.execute(sql)
    rows = target_cur.fetchall()
    for r in rows:

        # the the name of the current sample
        sql = "SELECT sample_name FROM samples WHERE pk_id=%s"
        target_cur.execute(sql, (r['fk_sample_id'], ))
        sample_name = target_cur.fetchone()[0]

        means = {}

        for lvl in ['t0', 't5', 't10', 't25', 't50', 't100', 't250']:

            # get all other members of this cluster - DO NOT INCLUDE IGNORED SAMPLES
            sql = "SELECT s.sample_name AS samplename FROM samples s, sample_clusters c WHERE s.ignore_sample=false AND s.ignore_zscore=false AND c.fk_sample_id=s.pk_id AND c."+lvl+"=%s AND c.fk_sample_id!=%s"
            target_cur.execute(sql, (r[lvl], r['fk_sample_id'], ))

            # if there are no other members
            if target_cur.rowcount <= 0:
                means[lvl] = None
            else:
                # get the distances to all other members and store in list
                other_dists = []
                rows2 = target_cur.fetchall()
                for sam in [r2['samplename'] for r2 in rows2]:
                    try:
                        other_dists.append(dm[sam][sample_name])
                    except KeyError:
                        try:
                            other_dists.append(dm[sample_name][sam])
                        except KeyError:
                            logging.warning("Distance from %s to %s is missing. Calculating ..", sample_name, sam)
                            missing_d = calculate_missing_distance(target_cur, sample_name, sam)
                            logging.warning("Distance is %i", missing_d)
                            other_dists.append(missing_d)
                            #raise SystemExit('Problem with distance matrix. No entry for %s and %s', s1, s2)
                # store mean of distances to other members in dict
                means[lvl] = sum(other_dists)/float(len(other_dists))

        # check we got all of them
        assert len(means.keys()) == 7

        # update table with mean of all distances of the current cluster to all other members of the sample's clusters
        sql = "UPDATE sample_clusters SET (t0_mean, t5_mean, t10_mean, t25_mean, t50_mean, t100_mean, t250_mean) = (%s, %s, %s, %s, %s, %s, %s) WHERE pk_id=%s"
        target_cur.execute(sql, (means['t0'], means['t5'], means['t10'], means['t25'], means['t50'], means['t100'], means['t250'], r['pk_id'], ))

    return 0

# --------------------------------------------------------------------------------------------------

def calculate_per_cluster_stats(target_cur, dm):
    """
    Calculate statistic for each cluster and add to cluster_stats table.

    Parameters
    ----------
    target_cur: obj
        database cursor
    dm: dict
        container for distance matrix

    Returns
    -------
    always 0
    """

    for lvl in ['t0', 't5', 't10', 't25', 't50', 't100', 't250']:

        # get all clusters on a level and the number of members WITHOUT IGNORE SAMPLES from database
        # also DON'T INCLUDE zscore_ignore SAMPLES
        sql = "SELECT c.%s, count(DISTINCT c.fk_sample_id) AS nof_members FROM sample_clusters c, samples s WHERE c.fk_sample_id = s.pk_id AND s.ignore_zscore=false AND s.ignore_sample=false GROUP BY c.%s" % (lvl, lvl)
        target_cur.execute(sql)
        rows = target_cur.fetchall()
        logging.info("Processing %i clusters at level %s", len(rows), lvl)
        for r in rows:

            # cluster 'name'
            clu = r[lvl]
            # number of members and number of pairwise dists
            nof_mems = r['nof_members']
            nof_pw_dists = ((nof_mems**2) - nof_mems) / 2

            # if we have only one member we have no pw dists and we cannot compute mean and ssd
            if nof_pw_dists > 0:
                # get all other members of the cluster DO NOT INCLUDE IGNORED SAMPLE IN COMPUTATION
                # also DON'T INCLUDE zscore_ignore SAMPLES
                sql = "SELECT s.sample_name FROM samples s, sample_clusters c WHERE s.ignore_sample=false AND s.ignore_zscore=false AND c.fk_sample_id=s.pk_id AND c." + lvl + "=%s"
                target_cur.execute(sql, (clu, ))
                clu_samples = [x['sample_name'] for x in target_cur.fetchall()]

                # check this complies with the count from above
                assert len(clu_samples) == nof_mems

                clu_dists = get_all_pw_dists(target_cur, clu_samples, dm)

                # check we get the expected number of distances back
                assert len(clu_dists) == nof_pw_dists

                mean_pw_dists = sum(clu_dists) / float(nof_pw_dists)
                ssd = sum([(x-mean_pw_dists)**2 for x in clu_dists])
                sd = sqrt(ssd / nof_pw_dists)
            else:
                mean_pw_dists = None
                sd = None

            sql = "INSERT INTO cluster_stats (cluster_level, cluster_name, nof_members, nof_pairwise_dists, mean_pwise_dist, stddev) VALUES (%s, %s, %s, %s, %s, %s)"
            target_cur.execute(sql, (lvl, clu, nof_mems, nof_pw_dists, mean_pw_dists, sd, ))

    return 0

# --------------------------------------------------------------------------------------------------

def get_all_pw_dists(target_cur, samples, dm):
    """
    Get all pairwise distances between the sampkes in the list.

    Parameters
    ----------
    samples: list
        list of sample names (MUST NOT contain ignored sample)
    dm: dict
        container for distance matrix

    Returns
    -------
    dists: list
        list of integers the distances
    raises: SystemExit
        if the distance matrix in the original database does not contain the distance
        between to samples that are both not ignored

    """

    dists = []
    for i, s1 in enumerate(samples):
        for j, s2 in enumerate(samples):
            if i>j:
                try:
                    dists.append(dm[s1][s2])
                except KeyError:
                    try:
                        dists.append(dm[s2][s1])
                    except KeyError:
                        logging.warning("Distance from %s to %s is missing. Calculating ..", s1, s2)
                        missing_d = calculate_missing_distance(target_cur, s1, s2)
                        logging.warning("Distance is %i", missing_d)
                        dists.append(missing_d)
                        #raise SystemExit('Problem with distance matrix. No entry for %s and %s', s1, s2)
    return dists

# --------------------------------------------------------------------------------------------------

def calculate_missing_distance(target_cur, s1, s2):

    sql = "SELECT pk_id FROM samples WHERE sample_name=%s"
    target_cur.execute(sql, (s1, ))
    s1_id = target_cur.fetchone()[0]

    sql = "SELECT pk_id FROM samples WHERE sample_name=%s"
    target_cur.execute(sql, (s2, ))
    s2_id = target_cur.fetchone()[0]

    # get list of contig ids from database
    sql = "SELECT pk_id FROM contigs"
    target_cur.execute(sql)
    rows = target_cur.fetchall()
    contig_ids = [r['pk_id'] for r in rows]

    d = {}
    for cid in contig_ids:
        t0 = time()
        target_cur.callproc("get_sample_distances_by_id", [s1_id, cid, [s2_id]])
        result = target_cur.fetchall()
        t1 = time()
        logging.info("Calculated %i distances on contig %i with 'get_sample_distances_by_id' in %.3f seconds",
                     len(result), cid, t1 - t0)

        # sum up if there are more than one contigs
        for res in result:
            if res[2] == None:
                res[2] = 0
            try:
                d[res[0]] += res[2]
            except KeyError:
                d[res[0]] = res[2]

    return d[s2_id]

# --------------------------------------------------------------------------------------------------

def read_distance_matrix(source_cur, dm):
    """
    Function to migrate the contig information.

    Parameters
    ----------
    source_cur: obj
        database cursor
    dm: dict
        container for distance matrix

    Returns
    -------
    no returns, but writes to dm
    """

    sql = "SELECT strain1, strain2, snp_dist FROM dist_matrix"
    source_cur.execute(sql)
    rows = source_cur.fetchall()
    for r in rows:
        try:
            dm[r['strain1']][r['strain2']] = r['snp_dist']
        except KeyError:
            dm[r['strain1']] = {r['strain2']: r['snp_dist']}

    return 0

# --------------------------------------------------------------------------------------------------

def migrate_contigs(source_cur, target_cur):
    """
    Function to migrate the contig information.

    Parameters
    ----------
    source_cur: obj
        database cursor
    target_cur: obj
        database cursor

    Returns
    -------
    all_contigs: set
        names of all contigs
    """

    all_contigs = set()

    for tbl in ['variants', 'ignored_pos']:
        sql = "SELECT DISTINCT contig FROM %s" % (tbl)
        source_cur.execute(sql)
        rows = source_cur.fetchall()
        all_contigs.update([row['contig'] for row in rows])

    for contig in all_contigs:
        sql = "INSERT INTO contigs (name) VALUES (%s)"
        target_cur.execute(sql, (contig, ))

    return all_contigs

# --------------------------------------------------------------------------------------------------

def migrate_variants(ref, source_cur, target_cur, all_sample_names):
    """
    Function to migrate the variant information.

    Parameters
    ----------
    ref: str
        name of reference in the database
    source_cur: obj
        database cursor
    target_cur: obj
        database cursor
    all_sample_names: set
        set with all samples names in db


    Returns
    -------
    no return,
    """

    logging.info("Reading variant and ign position information into memory.")

    # read contig variant and ignore position information into dicts
    dContigs = {}
    sql = "SELECT pk_id, name FROM contigs"
    target_cur.execute(sql)
    rows = target_cur.fetchall()
    dContigs = {r['name']: r['pk_id'] for r in rows}

    dVars = {}
    sql = "SELECT id, pos, var_base, contig FROM variants"
    source_cur.execute(sql)
    rows = source_cur.fetchall()
    for r in rows:
        dVars[r['id']] = {'pos': r['pos'] , 'var_base': r['var_base'], 'contig': r['contig']}

    dIgn = {}
    sql = "SELECT id, pos, contig FROM ignored_pos"
    source_cur.execute(sql)
    rows = source_cur.fetchall()
    for r in rows:
        dIgn[r['id']] = {'pos': r['pos'] , 'contig': r['contig']}

    # add reference to the database first and return the ignore positions as a set
    ref_ign_pos = add_ref_and_get_ign_pos(ref, source_cur, target_cur, dContigs, dIgn)
    if ref_ign_pos == None:
        logging.error("Could not get ref ign pos.")
        raise SystemExit("critical error")

    for sam in all_sample_names:

        # reference is already processed
        if sam == ref:
            continue

        logging.info("Reformatting variants for %s", sam)

        # get the sample id for this sample
        sql = "SELECT pk_id FROM samples WHERE sample_name=%s"
        target_cur.execute(sql, (sam, ))
        sample_id = target_cur.fetchone()[0]

        # get the infor from the variant table about this sample
        sql = "SELECT id, name, variants_id, ignored_pos FROM strains_snps WHERE name=%s ORDER BY id DESC"
        source_cur.execute(sql, (sam, ))
        rows = source_cur.fetchall()
        # if a samples has multiple sets of variants that are different to each other use the one with
        # the highest id, which is the first element in rows because of the above "ORDER BY id DESC"
        if len(rows) > 1:
            if len(rows) != rows.count(rows[0]):
                mess = "Multiple sets of different variant information found for sample %s. Using most recent." % (sam)
                logging.warning(mess)
        # let hope this is an ignore sample ...
        if len(rows) < 1:
            logging.warning("No variant information found for %s.", sam)
            continue

        # initialise a variant container for this sample with A, C, G, T
        sample_vars = {}
        for conname in dContigs.keys():
            sample_vars[conname] = {'A': set(), 'C': set(), 'G': set(), 'T': set()}

        r = rows[0]
        for vid in r['variants_id']:
            # get relavent infomation for this variant from the dict
            vb = dVars[vid]['var_base']
            p = dVars[vid]['pos']
            contig_name = dVars[vid]['contig']
            # add variant to set
            sample_vars[contig_name][vb].add(p)

        # create and add N pos set for each contig in one go
        for conname in dContigs.keys():
            sample_vars[conname]['N'] = set([dIgn[iid]['pos'] for iid in r['ignored_pos'] if dIgn[iid]['contig'] == conname])
            sample_vars[conname]['-'] = set()

        data = []
        for conname, contig_id in dContigs.iteritems():

            # remove ignore positions for the reference from all position sets in this sample
            for x in ['A', 'C', 'G', 'T', 'N', '-']:
                sample_vars[conname][x].difference_update(ref_ign_pos[conname])

            data.append((sample_id,
                         contig_id,
                         list(sample_vars[conname]['A']),
                         list(sample_vars[conname]['C']),
                         list(sample_vars[conname]['G']),
                         list(sample_vars[conname]['T']),
                         list(sample_vars[conname]['N']),
                         list(sample_vars[conname]['-'])))
        sql = "INSERT INTO variants (fk_sample_id, fk_contig_id, a_pos, c_pos, g_pos, t_pos, n_pos, gap_pos) VALUES (%s, %s, %s, %s, %s, %s, %s, %s)"
        target_cur.executemany(sql, data)

    logging.info("Added all variants to db.")

    return

# --------------------------------------------------------------------------------------------------

def add_ref_and_get_ign_pos(ref, source_cur, target_cur, dContigs, dIgn):
    """

    Parameters
    ----------
    ref: str
        reference name in db
    source_cur: obj
        database cursor
    target_cur: obj
        database cursor
    dContigs: dict
        this one: {r['name']: r['pk_id'] for r in rows}
    dIgn: dict
        dIgn[r['id']] = {'pos': r['pos'] , 'contig': r['contig']}
    Returns
    -------
    sample_vars: dict
        key: contig name, value: set of ign pos

    """

    logging.info("Reformatting variants for referemce %s", ref)

    # get the sample id for this sample
    sql = "SELECT pk_id FROM samples WHERE sample_name=%s"
    target_cur.execute(sql, (ref, ))
    sample_id = target_cur.fetchone()[0]

    # get the infor from the variant table about this sample
    sql = "SELECT id, ignored_pos FROM strains_snps WHERE name=%s ORDER BY id DESC"
    source_cur.execute(sql, (ref, ))
    rows = source_cur.fetchall()
    # if a samples has multiple sets of variants that are different to each other use the one with
    # the highest id, which is the first element in rows because of the above "ORDER BY id DESC"
    if len(rows) > 1:
        if len(rows) != rows.count(rows[0]):
            mess = "Multiple sets of different variant information found for sample %s. Using most recent." % (ref)
            logging.warning(mess)
    # let hope this is an ignore sample ...
    if len(rows) < 1:
        logging.error("No variant information found for %s.", ref)
        return None

    # initialise a variant container for this sample with A, C, G, T
    sample_vars = {}
    r = rows[0]

    # create and add N pos set for each contig in one go
    for conname in dContigs.keys():
        sample_vars[conname] = set([dIgn[iid]['pos'] for iid in r['ignored_pos'] if dIgn[iid]['contig'] == conname])

    data = []
    for conname, contig_id in dContigs.iteritems():
        data.append((sample_id,
                     contig_id,
                     list(),
                     list(),
                     list(),
                     list(),
                     list(sample_vars[conname]),
                     list()))
    sql = "INSERT INTO variants (fk_sample_id, fk_contig_id, a_pos, c_pos, g_pos, t_pos, n_pos, gap_pos) VALUES (%s, %s, %s, %s, %s, %s, %s, %s)"
    target_cur.executemany(sql, data)

    return sample_vars

# --------------------------------------------------------------------------------------------------

def migrate_samples_and_clusters(source_cur, target_cur):
    """
    Function to migrate the samples into the samples table of the new database.

    Parameters
    ----------
    source_cur: obj
        database cursor
    target_cur: obj
        database cursor

    Returns
    -------
    all_sample_names: set
        set with all samples names in db
    could raise SystemExit
    """

    # get set of all samples names as the union of what the names found in
    # strain_clusters, strain_stats, and strains_snps
    all_sample_names = set()

    for tbl in ['strain_clusters', 'strain_stats', 'strains_snps']:

        sql = "SELECT DISTINCT name FROM %s" % (tbl)
        source_cur.execute(sql)
        rows = source_cur.fetchall()
        all_sample_names.update([row['name'] for row in rows])

    logging.info("Found %i distinct sample names in source database.", len(all_sample_names))

    # for each of these names ...
    for sam in all_sample_names:

        # ... make an entry in the samples table and get the primary sample id
        sql = "INSERT INTO samples (sample_name) VALUES (%s) RETURNING pk_id"
        target_cur.execute(sql, (sam, ))
        sample_pkid = target_cur.fetchone()[0]

        # .. get the relevant information from strain stats
        sql = "SELECT time_of_upload, ignore, zscore_check FROM strain_stats WHERE name=%s ORDER BY time_of_upload ASC"
        source_cur.execute(sql, (sam, ))
        rows = source_cur.fetchall()

        # if there is more than one row for this sampke in strain_stats, print a warning
        if len(rows) > 1:
            logging.warning("More than one row found in strain_stats for sample %s. Only the most recent will be migrated.", sam)

        # update data in samples table -> overwritten by last entry in rows if more than one
        # last entry is most recent because of the "ORDER BY time_of_upload ASC" above
        ignore_flag = False
        for r in rows:
            ignore_flag = False
            zscore_flag = False
            if r['ignore'] != None:
                ignore_flag = True
            if r['zscore_check'] != None:
                zscore_flag = True
            sql = "UPDATE samples SET (date_added, ignore_sample, ignore_zscore) = (%s, %s, %s) WHERE pk_id=%s"
            target_cur.execute(sql, (r['time_of_upload'], ignore_flag, zscore_flag, sample_pkid, ))

        # this sample is ignored -> don't transfer clustering information
        if ignore_flag == True:
            continue

        # get clustering information from strain_clusters
        sql = "SELECT t0, t5, t10, t25, t50, t100, t250 FROM strain_clusters WHERE name=%s"
        source_cur.execute(sql, (sam, ))
        rows = source_cur.fetchall()
        # abandon all hope if there is conflicting information about the clustering
        if len(rows) > 1:
            if len(rows) != rows.count(rows[0]):
                mess = "Multiple sets of different clustering information found for sample %s. Could not migrate database." % (sam)
                logging.error(mess)
                raise SystemExit(mess)
        # let hope this is an ignore sample ...
        if len(rows) < 1:
            logging.warning("No clustering information found for %s.", sam)
            continue

        # else enter clustering information in sample_clusters
        r = rows[0]
        sql = "INSERT INTO sample_clusters (fk_sample_id, t0, t5, t10, t25, t50, t100, t250) VALUES (%s, %s, %s, %s, %s, %s, %s, %s)"
        target_cur.execute(sql, (sample_pkid, r['t0'], r['t5'], r['t10'], r['t25'], r['t50'], r['t100'], r['t250'], ))

    return all_sample_names

# --------------------------------------------------------------------------------------------------

if __name__ == '__main__':
    sys.exit(main())
