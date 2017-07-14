#!/usr/bin/env python
"""

"""
import sys
import argparse
import math

__version__= '0.1'
__date__= '20Apr2017'
__author__ = 'ulf.schaefer@phe.gov.uk'

# --------------------------------------------------------------------------------------------------

class ClusterStatsError(Exception):
    pass

# --------------------------------------------------------------------------------------------------

class ClusterStats(object):
    '''

    '''
    def __init__(self, **kwargs):
        """
        Constructor of the class
        Required parameter:
            - number of members (int)
        Accepts either:
            - pairwise distances (list of ints)
            OR
            - stddev of pw dists (float)
            - mean of pwdists(float)
        """

        self.members = int(kwargs['members'])
        self.nof_pw_dists = None
        self.mean_pw_dist = None
        self.stddev_pw_dist = None
        self.variance_pw_dist = None

        if kwargs.has_key('dists'):

            self.nof_pw_dists = len(kwargs['dists'])

            if self.nof_pw_dists != (self.members * (self.members-1))/2.0:
                raise ClusterStatsError("Nof members and nof distances inconsistent.")

            if self.nof_pw_dists > 0:
                self.mean_pw_dist = float(sum(kwargs['dists'])/self.nof_pw_dists)
                #calculate variance and stddev
                x = []
                for d in kwargs['dists']:
                    x.append((d - self.mean_pw_dist)**2.0)

                self.variance_pw_dist = sum(x)/len(x)
                self.stddev_pw_dist = math.sqrt(self.variance_pw_dist)
            else:
                # this is for one member only clusters
                self.mean_pw_dist = None
                self.variance_pw_dist = None
                self.stddev_pw_dist = None

        elif kwargs.has_key('stddev') and kwargs.has_key('mean'):
            self.nof_pw_dists = (self.members * (self.members-1))/2.0
            self.mean_pw_dist = float(kwargs['mean'])
            self.stddev_pw_dist = float(kwargs['stddev'])
            self.variance_pw_dist = self.stddev_pw_dist**2.0
        else:
            raise ClusterStatsError("kwargs combination passed to constructor not valid.")

    # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

    def __str__(self):
        '''
        overwrite __str__() to tell it how to print itself
        '''
        return "members: %s, nof_pw_dists: %s, mean_pw_dist: %s, stddev_pw_dist: %s, variance_pw_dist: %s" \
                %(self.members,
                  self.nof_pw_dists,
                  self.mean_pw_dist,
                  self.stddev_pw_dist,
                  self.variance_pw_dist)

    # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

    def add_member(self, new_dists):
        '''
        add a new member to the object
        '''

        assert len(new_dists) == self.members

        if self.members > 1:
            for nd in new_dists:
                # remember the mean before updating
                prev_m = self.mean_pw_dist

                # proxy for the sum of all pw dists
                sm = self.mean_pw_dist * float(self.nof_pw_dists)
                # add one new dist to the sum
                new_sum = sm + nd
                # we have one more dist at this point
                self.nof_pw_dists += 1
                # new mean is new_sum over new nof dists
                self.mean_pw_dist = new_sum / self.nof_pw_dists

                # update the variance
                # see https://math.stackexchange.com/questions/775391
                N = self.nof_pw_dists
                a = (N - 1) * (self.variance_pw_dist)
                b = (nd - self.mean_pw_dist) * (nd - prev_m)
                self.variance_pw_dist = (a + b) / N
                self.stddev_pw_dist = math.sqrt(self.variance_pw_dist)
        else:
            self.nof_pw_dists = 1
            self.mean_pw_dist = float(new_dists[0])
            self.stddev_pw_dist = 0.0
            self.variance_pw_dist = 0.0

        self.members += 1

# --------------------------------------------------------------------------------------------------
