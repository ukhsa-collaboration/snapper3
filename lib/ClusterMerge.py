#!/usr/bin/env python
"""

"""
import sys
import argparse
import math

__version__= '0.1'
__date__= '14Jul2017'
__author__ = 'ulf.schaefer@phe.gov.uk'

# --------------------------------------------------------------------------------------------------

class ClusterMerge(object):
    '''
    Class to help with cluster merging.

    '''
    def __init__(self, **kwargs):
        """
        Constructor of the class. Requires two kwarg:
            - level (int)
            - clusters (list of int)
        """

        self.level = int(kwargs['level'])
        self.t_level = 't%i' % self.level
        self.org_clusters = kwargs['clusters']
        self.final_name = None
        self.final_members = None
        self.stats = None

    # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

    def __str__(self):
        '''
        overwrite __str__() to tell it how to print itself
        '''
        return "level: %s, clusters: %s" \
                %(self.t_level,
                  self.org_clusters)

    # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

# --------------------------------------------------------------------------------------------------

