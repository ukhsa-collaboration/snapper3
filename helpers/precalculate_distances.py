import sys
import os
import argparse
import logging
import json
import gzip
import psycopg2
from psycopg2.extras import DictCursor

from operator import itemgetter

__version__ = '0.1'
__date__ = '05Sep2018'
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

    parser.add_argument("--connstring",
                        "-c",
                        type=str,
                        metavar="CONNECTION",
                        required=True,
                        dest="db",
                        help="Connection string for db. REQUIRED.")

    parser.add_argument("--jsonvariants",
                        "-j",
                        type=str,
                        metavar="JSONFILE",
                        required=True,
                        dest="vars",
                        help="Variants for this sample as PHEnix make jsons file. REQUIRED.")

    parser.add_argument("--sample_name",
                        "-s",
                        type=str,
                        metavar="SAMPLENAME",
                        default=None,
                        dest="sam_name_in",
                        help="Sample name. Default: Name of json file before the 1st dot.")

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

    logging.basicConfig(format="[%(asctime)s] %(levelname)s: %(message)s", level=logging.DEBUG)

    if oArgs.sam_name_in == None:
        oArgs.sam_name_in = os.path.basename(oArgs.vars).split('.')[0]

    open_func = gzip.open if oArgs.vars.endswith('.gz') == True else open

    data = None
    with open_func(oArgs.vars) as f:
        data = json.load(f)

    # print data['positions'][u'gi|206707319|emb|AM933172.1|'].keys()
    for cnt in data['positions'].keys():
        for nt in data['positions'][cnt].keys():
            data['positions'][cnt][nt] = set(data['positions'][cnt][nt])

    try:
        # open source db
        conn = psycopg2.connect(oArgs.db)
        cur = conn.cursor(cursor_factory=psycopg2.extras.DictCursor)

        sql = "SELECT pk_id, name FROM contigs"
        cur.execute(sql)
        rows = cur.fetchall()
        contigs = {r['pk_id']: r['name'] for r in rows}

        logging.debug("%i contigs found.", len(contigs))

        # get all samples that are NOT ignored and that have a SNP address
        sql = "SELECT c.fk_sample_id, s.sample_name FROM samples s, sample_clusters c WHERE s.pk_id=c.fk_sample_id AND s.ignore_sample IS FALSE"
        cur.execute(sql)
        rows = cur.fetchall()
        samples = {r['fk_sample_id']: r['sample_name'] for r in rows}

        logging.debug("%i samples found.", len(samples))

        sql = "SELECT fk_sample_id, fk_contig_id, a_pos, c_pos, g_pos, t_pos, n_pos, gap_pos FROM variants WHERE fk_sample_id IN %s"
        cur.execute(sql, (tuple(samples.keys()), ))
        variants = {}
        for row in cur.fetchall():
            try:
                _ = variants[row['fk_sample_id']]
            except KeyError:
                variants[row['fk_sample_id']] = {}
            variants[row['fk_sample_id']][row['fk_contig_id']] = {'A': set(row['a_pos']),
                                                                  'C': set(row['c_pos']),
                                                                  'G': set(row['g_pos']),
                                                                  'T': set(row['t_pos']),
                                                                  'N': set(row['n_pos'] + row['gap_pos'])}

        logging.debug("%i variant sets found.", len(variants))

    except psycopg2.Error as e:
        logging.error("Database reported error: %s" % (str(e)))
    finally:
        # close all dbs
        cur.close()
        conn.close()

    logging.debug("Calculating distances...")

    distances = []
    for samid, sam_name in samples.iteritems():
        d = 0
        for c_id, c_nme in contigs.iteritems():
            d += len(((variants[samid][c_id]['A'] ^ data['positions'][c_nme]['A']) |
                      (variants[samid][c_id]['C'] ^ data['positions'][c_nme]['C']) |
                      (variants[samid][c_id]['G'] ^ data['positions'][c_nme]['G']) |
                      (variants[samid][c_id]['T'] ^ data['positions'][c_nme]['T'])) \
                     - (variants[samid][c_id]['N'] | data['positions'][c_nme]['N'] | data['positions'][c_nme]['-']))

        distances.append((samid, sam_name, d))

    distances.sort(key=lambda x: x[2])

    resu = {'sample_name': oArgs.sam_name_in, 'distances': distances}
    outfilename = '%s.distances.json.gz' % (oArgs.sam_name_in)
    with gzip.open(outfilename, 'w') as outfile:
        json.dump(resu, outfile)

    logging.debug("Written %i distances to file %s", len(distances), outfilename)

    return 0

# end of main --------------------------------------------------------------------------------------

if __name__ == '__main__':
    sys.exit(main())
