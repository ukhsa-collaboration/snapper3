"""
Module for accessing the SnapperDB3 database.

:author: ulf.schaefer@phe.gov.uk
:created: 25 July 2016
"""

import logging

import psycopg2
from psycopg2.extras import DictCursor


from lib.distances import get_distances, get_relevant_distances
from lib.utils import get_closest_threshold

# --------------------------------------------------------------------------------------------------

class SnapperDBInterrogationError(Exception):
    pass

# --------------------------------------------------------------------------------------------------

class SnapperDBInterrogation(object):
    """
    Class to wrap a database connection to a Snapper3 db and expose
    the required functionality.

    """

    # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

    def __init__(self, **kwargs):
        """
        Constructor.

        Parameters:
        -----------

        """

        self.connstring = None

        if kwargs.has_key('conn_string'):
            self.connstring = kwargs['conn_string']
        elif all(kwargs.has_key(x) for x in ["host", "dbname", "user", "password"]) == True:
            self.connstring = "host='%s' dbname='%s' user='%s' password='%s'" % \
                              (kwargs["host"],
                               kwargs["dbname"],
                               kwargs["user"],
                               kwargs["password"])
        else:
            raise SnapperDBInterrogationError("kwargs combination passed to constructor not valid.")

    # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

    def _connect(self):
        """
        **PRIVATE**

        Connect to the external database for querying.
        """

        self.conn = psycopg2.connect(self.connstring)
        self.cur = self.conn.cursor(cursor_factory=DictCursor)

    # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

    def _close(self):
        """
        **PRIVATE**

        Close the connection to the external database.

        There is no "commit" here because we will not be writing anything to the database.

        """

        if not self.cur.closed:
            self.cur.close()
        if self.conn.closed == 0:
            self.conn.close()

    # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

    def __enter__(self):
        """
        Enter function to use with context manager.
        """

        # open connection for use with context managers.
        if not hasattr(self, "cur") or self.cur.closed:
            self._connect()

        return self

    # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

    def __exit__(self, exit_type, value, traceback):
        """
        Exit function to use with context manager.
        """

        # Close the connection for use with context managers.
        self._close()

        return False

    # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

    def get_closest_samples(self, sam_name, neighbours, levels=[0, 5, 10, 25, 50, 100, 250]):
        """
        Get the closest n samples.

        Parameters
        ----------
        sam_name: str
            name of the query sample
        neighbours: int
            number on neighbours
        levels: list of int
            default: [0, 5, 10, 25, 50, 100, 250]
            better don't change it

        Returns
        -------
        result_samples: list of tuples
            sorted [(sample_name, distance), (sample_name, distance), (sample_name, distance), ..]
        """

        # get the snp address of the query sample
        sql = "SELECT s.pk_id, c.t0, c.t5, c.t10, c.t25, c.t50, c.t100, c.t250 FROM sample_clusters c, samples s WHERE s.pk_id=c.fk_sample_id AND s.sample_name=%s"
        self.cur.execute(sql,(sam_name, ))
        if self.cur.rowcount < 1:
            raise SnapperDBInterrogationError("No clustering information found for sample %s" % (sam_name))
        row = self.cur.fetchone()
        snad = [row['t0'], row['t5'], row['t10'], row['t25'], row['t50'], row['t100'], row['t250']]
        samid = row['pk_id']

        close_samples = set()
        id2name = {}
        for clu, lvl in zip(snad, levels):

            t_lvl = 't%i' % (lvl)
            sql = "SELECT s.sample_name, c.fk_sample_id FROM sample_clusters c, samples s WHERE c."+t_lvl+"=%s AND s.pk_id=c.fk_sample_id"
            self.cur.execute(sql, (clu, ))
            rows = self.cur.fetchall()
            for r in rows:
                id2name[r['fk_sample_id']] = r['sample_name']
                if r['fk_sample_id'] != samid:
                    close_samples.add(r['fk_sample_id'])

            logging.info("Number of samples in same %s cluster: %i.", t_lvl, len(close_samples))

            if len(close_samples) >= neighbours:
                break

        distances = None
        if len(close_samples) < neighbours:
            distances = get_relevant_distances(self.cur, samid)
            sql = "SELECT pk_id, sample_name FROM samples"
            self.cur.execute(sql)
            id2name = {r['pk_id']: r['sample_name'] for r in self.cur.fetchall()}
        else:
            distances = get_distances(self.cur, samid, list(close_samples))
        result_samples = distances[:neighbours]

        for (sa, di) in distances[neighbours:]:
            if di == result_samples[-1][1]:
                result_samples.append((sa, di))

        result_samples = [(id2name[sa], di) for (sa, di) in result_samples]
        return result_samples

    # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

    def get_samples_below_threshold(self, sam_name, dis, levels=[0, 5, 10, 25, 50, 100, 250]):
        """
        Get all samples that are below or equal to a given distance from the query sample.

        Parameters
        ----------
        sam_name: str
            query sample name
        dis: int
            distance threshold
        levels: list of ints
            default: [0, 5, 10, 25, 50, 100, 250]
            better don't change it

        Returns
        -------
        result_samples: list of tuples
            sorted [(sample_name, distance), (sample_name, distance), (sample_name, distance), ..]
        """

        # get the snp address of the query sample
        sql = "SELECT s.pk_id, c.t0, c.t5, c.t10, c.t25, c.t50, c.t100, c.t250 FROM sample_clusters c, samples s WHERE s.pk_id=c.fk_sample_id AND s.sample_name=%s"
        self.cur.execute(sql,(sam_name, ))
        if self.cur.rowcount < 1:
            raise SnapperDBInterrogationError("No clustering information found for sample %s" % (sam_name))
        row = self.cur.fetchone()
        snad = [row['t0'], row['t5'], row['t10'], row['t25'], row['t50'], row['t100'], row['t250']]
        samid = row['pk_id']

        ct = get_closest_threshold(dis)
        t_ct = 't%i' % (ct)
        cluster = snad[levels.index(ct)]

        id2name = {}
        sql = "SELECT s.sample_name AS samname, c.fk_sample_id AS samid FROM sample_clusters c, samples s WHERE c."+t_ct+"=%s AND s.pk_id=c.fk_sample_id"
        self.cur.execute(sql, (cluster, ))
        rows = self.cur.fetchall()
        neighbours = []
        for r in rows:
            id2name[r['samid']] = r['samname']
            if r['samid'] != samid:
                neighbours.append(r['samid'])

        if len(neighbours) <= 0:
            logging.info("No samples found this close to the query sample.")
            return []
        else:
            logging.info("Calculating distances to %i samples in the same %s cluster %s.", len(neighbours), t_ct, cluster)
            distances = get_distances(self.cur, samid, neighbours)
            result_samples = [(id2name[s], d) for (s, d) in distances if d <= dis]
            return result_samples

    # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

    def get_snp_address(self, sam_name):
        """
        Get the snp address for the sample.

        Parameters
        ----------
        sam_name: str
            name of the query sample

        Returns
        -------
        snad: str
            "1-2-3-4-5-6-7" if successful else None
        """

        sql = "SELECT c.t0, c.t5, c.t10, c.t25, c.t50, c.t100, c.t250 FROM sample_clusters c, samples s WHERE s.pk_id=c.fk_sample_id AND s.sample_name=%s"
        self.cur.execute(sql, (sam_name, ))
        if self.cur.rowcount < 1:
            raise SnapperDBInterrogationError("No clustering information found for sample %s" % (sam_name))
        elif self.cur.rowcount > 1:
            raise SnapperDBInterrogationError("Too much clustering information found for sample %s" % (sam_name))
        else:
            row = self.cur.fetchone()
            return "%i-%i-%i-%i-%i-%i-%i" % (row['t250'], row['t100'], row['t50'], row['t25'], row['t10'], row['t5'], row['t0'])

# --------------------------------------------------------------------------------------------------
