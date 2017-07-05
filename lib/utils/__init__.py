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

def get_closest_threshold(i, levels=[0, 5, 10, 25, 50, 100, 250]):
    """

    Parameters
    ----------

    Returns
    -------
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
