import logging

import argparse
from argparse import RawTextHelpFormatter
import psycopg2
from psycopg2.extras import DictCursor

from lib.utils import read_fasta
import lib.alignment as align

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

    return r'''For a given set of sample names get either a SNP only or a whole genome
alignment. There are various filtering options at your disposal, which are identical to
Phenix's vcf2fasta.'''

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

    def positive_float(value):
        """Make type for float 0<x<1."""
        x = float(value)
        if not 0.0 <= x <= 1.0:
            raise argparse.ArgumentTypeError("%r not in range [0.0, 1.0]" % x)
        return x

    args = argparse.ArgumentParser(description=get_desc(), formatter_class=RawTextHelpFormatter)

    args.add_argument("--connstring",
                      "-c",
                      type=str,
                      metavar="CONNECTION",
                      required=True,
                      dest="db",
                      help="Connection string for db. REQUIRED")

    args.add_argument("--samples",
                      "-s",
                      type=str,
                      required=True,
                      nargs='+',
                      help="List of sample names to process. Separated by blank.")

    args.add_argument("--out",
                      "-o",
                      required=True,
                      help="Path to the output FASTA file.")

    args.add_argument("--reference",
                      type=str,
                      required=True,
                      help="Path to reference specified (FASTA). REQUIRED.")

    args.add_argument("--name-of-ref-in-db",
                      type=str,
                      required=True,
                      help="The name of the reference in the database. REQUIRED.")

    args.add_argument("--column-Ns",
                      type=positive_float,
                      help="Keeps columns with fraction of Ns below specified threshold.")
    args.add_argument("--column-gaps",
                      type=positive_float,
                      help="Keeps columns with fraction of Ns below specified threshold.")

    args.add_argument("--sample-Ns",
                      help="Keeps samples with fraction of Ns below specified threshold or put 'auto'." + \
                           "Fraction expressed as fraction of genome. Requires --reflength or --reference.")
    args.add_argument("--sample-gaps",
                      help="Keeps samples with fraction of gaps below specified threshold or put 'auto'." + \
                           "Fraction expressed as fraction of genome. Requires --reflength or --reference.")
    args.add_argument("--sample-Ns-gaps-auto-factor",
                      default=2.0,
                      type=float,
                      help="When using 'auto' option for --sample-gaps or --sample-Ns, remove sample that have" + \
                           "gaps or Ns this many times above the stddev of all samples. [Default: 2.0]")
    args.add_argument("--snp-address",
                      action='store_true',
                      help="Annotate fasta sample header with SNP address where available. [Default: don't]")

    group_b = args.add_mutually_exclusive_group()
    group_b.add_argument("--remove-invariant-npos",
                         action='store_true',
                         help="Remove all positions that are invariant apart from N positions.")
    group_b.add_argument("--whole-genome",
                         action='store_true',
                         help="Whole genome will be in the alignment. [Default: Variant positions only.]")

    group_c = args.add_mutually_exclusive_group()
    group_c.add_argument("--include",
                         help="Only include positions in BED file in the FASTA")
    group_c.add_argument("--exclude",
                         help="Exclude any positions specified in the BED file.")

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

    logging.debug("Args received: %s", str(args))

    # do some additional args checking
    if args['sample_Ns']:
        if args['sample_Ns'] != 'auto':
            try:
                args['sample_Ns'] = float(args['sample_Ns'])
                if not 0.0 <= args['sample_Ns'] <= 1.0:
                    raise TypeError
            except TypeError:
                logging.error("Please put either 'auto' or a float [0.0, 1.0] in sample-Ns.")
                return 1

    if args['sample_gaps']:
        if args['sample_gaps'] != 'auto':
            try:
                args['sample_gaps'] = float(args['sample_gaps'])
                if not 0.0 <= args['sample_gaps'] <= 1.0:
                    raise TypeError
            except TypeError:
                logging.error("Please put either 'auto' or a float [0.0, 1.0] in sample-gaps.")
                return 1

    ref_seq = {}
    with open(args["reference"]) as fp:
        ref_seq = read_fasta(fp)

    args["reference"] = ref_seq
    args['reflength'] = sum([len(x) for x in ref_seq.values()])

    nof_samples = len(args['samples'])
    logging.info("Processing %i samples plus the reference.", nof_samples)

    all_contig_data = align.get_data_from_db(args['db'], args['samples'], args['name_of_ref_in_db'])
    if all_contig_data == None:
        logging.error("There was a problem getting the data from the database.")
        return 1

    dSnads = None
    if args['snp_address'] == True:
        dSnads = align.get_snp_addresses(args['db'], args['samples'])
        if dSnads == None:
            logging.error("There was a problem getting the snp addresses from the database.")
            return 1

    if align.add_reference_data(args['reference'], all_contig_data) == None:
        logging.error("There was a problem adding reference data to the data structure.")
        return 1

    # check that there is no conflicting ref bases by verifying that the
    # intersection between two ref bases sets of positions is always empty
    for (contig, data) in all_contig_data.iteritems():
        ref_bases = data['reference'].keys()
        for i in range(0, len(ref_bases)):
            for j in range(0, len(ref_bases)):
                if i < j:
                    x = ref_bases[i]
                    y = ref_bases[j]
                    itscn = data['reference'][x] & data['reference'][y]
                    assert len(itscn) == 0, "FATAL ERROR: Two different ref bases for the same position: %s on contig %s" % (str(itscn), contig)

    """
    all_contig_data looks like this:
    {'gi|194097589|ref|NC_011035.1|': {'211696_H14256028001': {'-': set([12673,
                                                                         12674,
                                                                         12675,
                                                                         12676,
                                                                         12677]),
                                                                'A': set([52677,
                                                                          52678,
                                                                          52679,
                                                                          52680])
                                                                'C': set([52677,
                                                                          52678,
                                                                          52679,
                                                                          52680])
                                                                 ...
                                                               }
                                       }
                                       {'211697_H14333548001': {'-': set([12673,
                                                                         12674,
                                                                         12675,
                                                                         12676,
                                                                         12677]),
                                                                'A': set([52677,
                                                                          52678,
                                                                          52679,
                                                                          52680])
                                                                'C': set([52677,
                                                                          52678,
                                                                          52679,
                                                                          52680])
                                                                 ...
                                        ...

                                                               }
        ...
    }
    """

    # insert any filtering here ...
    logging.info("Filtering alignment.")

    if args["exclude"] or args["include"]:
        align.process_bed_file(args, all_contig_data)

    if args['remove_invariant_npos'] == True:
        logging.info("Removing invariant N positions.")
        for (contig, data) in all_contig_data.iteritems():
            # get all positions that are N and all others in all samples except the reference
            n_pos = set()
            var_pos = set()
            for sam in data.keys():
                if sam != 'reference':
                    for nuc in data[sam].keys():
                        if nuc == 'N':
                            n_pos.update(data[sam][nuc])
                        else:
                            var_pos.update(data[sam][nuc])
            # get positions that are onlys in nothing else in any sample
            only_n_pos = n_pos.difference(var_pos)
            # remove those positions from all samples including the reference
            for sam in data.keys():
                for nuc in data[sam].keys():
                    data[sam][nuc].difference_update(only_n_pos)

    if args['sample_Ns']:
        align.remove_samples(args, 'sample_Ns', 'N', all_contig_data)

    if args['sample_gaps']:
        align.remove_samples(args, 'sample_gaps', '-', all_contig_data)

    if args['column_Ns']:
        align.remove_columns(args['column_Ns'], 'N', all_contig_data)

    if args['column_gaps']:
        align.remove_columns(args['column_gaps'], '-', all_contig_data)

    # finished filtering
    logging.info("Filtering complete.")

    # output per sample stats
    logging.info("Writing per sample stats.")
    align.output_per_sample_stats(all_contig_data)

    # output now
    dSeqs = {}
    for (contig, data) in all_contig_data.iteritems():
        dAlign = {}
        # get all positions
        if args["whole_genome"]:
            # all positions for whole contig when reference is required
            all_pos = {i + 1: i for i in range(len(args["reference"][contig]))}
        else:
            # get all positions for variant positions only else
            all_pos = set()
            for nuc in data['reference'].keys():
                all_pos.update(data['reference'][nuc])
            # all_pos is now a list ... and then a dict
            all_pos = sorted(all_pos)
            all_pos = {all_pos[i]: i for i in range(len(all_pos))}

        # 'initialise' sequence
        for sample_name in data.keys():
            if args["whole_genome"]:
                # this is effectively a deepcopy, otherwise we're writing on the same
                # copy of the reference string for all samples
                # dAlign[sample_name] = args["reference"][contig]
                dAlign[sample_name] = []
                for s in args["reference"][contig]:
                    dAlign[sample_name].append(s)
            else:
                # initialies with 0's
                dAlign[sample_name] = ['0'] * len(all_pos)
                # set all bases to reference
                for nuc in data['reference'].keys():
                    for i in data['reference'][nuc]:
                        seq_pos = all_pos[i]
                        dAlign[sample_name][seq_pos] = nuc

            # overwrite reference positions where necessary
            for nuc in data[sample_name].keys():
                for i in data[sample_name][nuc]:
                    seq_pos = all_pos[i]
                    dAlign[sample_name][seq_pos] = nuc

            seq = ''.join(dAlign[sample_name])
            try:
                dSeqs[sample_name] += seq
            except KeyError:
                dSeqs[sample_name] = seq

    # write to file
    with open(args["out"], "w") as fp:
        # write seqs to file
        for name, seq in dSeqs.iteritems():
            if dSnads != None and dSnads.has_key(name):
                fp.write(">%s %s\n%s\n" % (name, dSnads[name], seq))
            else:
                fp.write(">%s\n%s\n" % (name, seq))

    return 0

# --------------------------------------------------------------------------------------------------

if __name__ == "__main__":
    exit(main(vars(get_args().parse_args())))
