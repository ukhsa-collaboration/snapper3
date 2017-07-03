"""
File contains some utilities that are used in snapper v3.

author: ulf.schaefer@phe.gov.uk

"""

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
