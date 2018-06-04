"""
File contains some utilities that are used in snapper v3.

author: ulf.schaefer@phe.gov.uk

"""

import gzip
import logging
import json

# --------------------------------------------------------------------------------------------------

def check_json_format(data):
    """
    Checks that one of the top level keys in the dict is called 'positions'
    and that the dict in there contains all nucleotides for all contigs.

    Parameters
    ----------
    data: dict
        input data read from json file

    Returns
    -------
    op: boolean
        the result of the check (True=fine, False=not fine)
    """

    op = False
    if data.has_key('positions') == False:
        return op

    contigs = data['positions'].keys()
    op = all([sorted(data['positions'][c].keys()) == [u'-', u'A', u'C', u'G', u'N', u'T'] for c in contigs])

    return op

# --------------------------------------------------------------------------------------------------

def get_closest_threshold(i, levels=[0, 5, 10, 25, 50, 100, 250]):
    """
    Get the closest snp address threshold level for a given distance.

    Parameters
    ----------
    i: int
        distance

    Returns
    -------
    levels[x]: int
        closest threshold

    """

    x = 0
    try:
        while levels[x] < i:
            x += 1
    except IndexError:
        # this happens when the closest sample is >250 away
        return None
    return levels[x]

# --------------------------------------------------------------------------------------------------

def get_the_data_from_the_input(args):
    """

    Parameters
    ----------

    Returns
    -------
    """

    open_func = gzip.open if args['input'].endswith('.gz') == True else open

    data = None
    if args['format'] == 'json':
        try:
            with open_func(args['input']) as data_file:
                try:
                    data = json.load(data_file)
                except ValueError:
                    logging.error("Data in %s is corrupted.", args['input'])
                    return None
        except IOError:
            logging.error("Could not open file %s", args['input'])
            return None
        if check_json_format(data) == False:
            logging.error("Data in %s is not in the correct format. Please use the latest version of Phenix to make this file from a filtered vcf.", args['input'])
            return None
    elif args['format'] == 'fasta':
        dRef = None
        dInp = None
        try:
            with open_func(args['input']) as fasta:
                dInp = read_fasta(fasta)
        except IOError:
            logging.error("Could not open file %s", args['input'])
            return None

        try:
            open_func = gzip.open if args['reference'].endswith('.gz') == True else open
            with open_func(args['reference']) as ref:
                dRef = read_fasta(ref)
        except IOError:
            logging.error("Could not open file %s", args['reference'])
            return None

        data = get_data_from_seqs(dInp, dRef)
        if data == None:
            logging.error("There was a problem with the sequences.")
            return None

    else:
        logging.error("Unknown data format: %s", args['format'])
        return None

    return data

# --------------------------------------------------------------------------------------------------

def get_data_from_seqs(sam, ref):
    """
    Get a dictionary from the fasta sequence that looks the same as the one in the Phenix-made json.

    Parameters
    ----------
    sam: dict
        sequence dict from fasta for sample
    ref: dict
        sequence dict from fasta for reference


    Returns
    -------
    data: dict
        this is the data, innit?
    None if fails
    """

    data = {'positions': {}}

    sample_contigs = sorted(sam.keys())
    refere_contigs = sorted(ref.keys())

    if len(sample_contigs) != len(refere_contigs):
        logging.error("Not the same number of contigs in sample and ref. Is this the right reference?")
        return None

    # if there is more than one contig, the names need to match
    if len(refere_contigs) > 1:
        if refere_contigs != sample_contigs:
            logging.error("Contigs not the same between sample and ref. Is this the right reference?")
            return None
        if [len(ref[x]) for x in refere_contigs] != [len(sam[x]) for x in refere_contigs]:
            logging.error("Contig lengths not the same between sample and ref. Is this the right reference?")
            return None
    else:
        # else only the length needs to match
        reflen = len(ref[refere_contigs[0]])
        samlen = len(sam[sample_contigs[0]])
        if samlen != reflen:
            logging.error("Sequence length not the same between sample and ref. Is this the right reference?")
            return None

    for (r_con, s_con) in zip(refere_contigs, sample_contigs):
        data['positions'][r_con] = {u'-': set(), u'A': set(), u'C': set(), u'G': set(), u'N': set(), u'T': set()}
        ref_seq = ref[r_con].upper()
        sam_seq = sam[s_con].upper()

        for i, (r, s) in enumerate(zip(ref_seq, sam_seq)):
            if r != s:
                try:
                    data['positions'][r_con][s].add(i+1)
                except KeyError:
                    logging.error("Unknown character in sample sequence. Only [A,C,G,T,N,-] are allowed.")
                    return None

    return data

# --------------------------------------------------------------------------------------------------

def read_fasta(f):
    """
    Read fasta into a dict. Use header line without > up to 1st ' ' as keys.

    Parameters
    ----------
    f: file handle
        fasta file

    Returns
    -------
    d: dict
        {header: seq, ...}
    """

    d = {}
    sSeq = ""
    sHeader = ""
    for sLine in f:
        sLine = sLine.strip()
        if sLine.startswith(">"):
            if len(sSeq) > 0:
                d[sHeader] = sSeq
                sSeq = ""
            sHeader = sLine[1:].split(" ")[0]
            continue
        sSeq = sSeq + sLine
    d[sHeader] = sSeq
    return d

# --------------------------------------------------------------------------------------------------

def get_all_cluster_members(cur, c, t):
    """
    Get all member of a given cluster.

    Parameters
    ----------
    cur: obj
        database cursor
    c: int
        cluster name
    t: str
        cluster threshold t0 or t5 or t10 or t25 ...

    Returns
    -------
    neighbours: list
        list of sample ids
    """

    sql = "SELECT fk_sample_id AS samid FROM sample_clusters WHERE "+t+"=%s"
    cur.execute(sql, (c, ))
    rows = cur.fetchall()
    neighbours = []
    for r in rows:
        neighbours.append(r['samid'])
    return neighbours

# --------------------------------------------------------------------------------------------------

def calculate_nless_n50(data, fasta):
    """
    Calculate the n50 of N-less sequence.

    Parameters
    ----------
    data: dict
        as returned from get_data_from_seqs()
    fasta: str
        fasta input file

    Returns
    -------
    """

    with open(fasta) as f:
        dInp = read_fasta(f)

    lens = []
    for con, condata in data['positions'].iteritems():
        n_list = list(condata['N'])
        n_list.sort()

        # if there are no Ns
        if len(n_list) <= 0:
            lens.append(len(dInp[con]))
        else:
            for i in range(len(n_list)):
                if i == 0:
                    lens.append(n_list[i] - 1)
                if i == len(n_list) - 1:
                    lens.append(len(dInp[con]) - n_list[i])
                else:
                    l = (n_list[i] - n_list[i-1]) - 1
                    lens.append(l)

    lens = [x for x in lens if x > 0]
    if len(lens) <= 0:
        n50 = 0
    else:
        lens.sort(reverse=True)
        sm = sum(lens)
        x = 0
        i = 0
        while x < sm*0.5:
            x += lens[i]
            i += 1
        n50 = lens[i-1]

    return n50

# --------------------------------------------------------------------------------------------------
