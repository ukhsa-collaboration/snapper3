"""
Module for accessing the SnapperDB3 database.

:author: ulf.schaefer@phe.gov.uk
:created: 25 July 2016
"""

import logging
import subprocess
import tempfile
import os
import shutil

import psycopg2
from psycopg2.extras import DictCursor

from lib.distances import get_distances, get_relevant_distances, get_distance_matrix
from lib.distances import get_distances_fusion, get_relevant_distances_fusion
from lib.distances import get_distance_matrix_fusion
from lib.utils import get_closest_threshold

import get_alignment

HAVE_BIOPYTHON = True
try:
    from Bio import Phylo
    from Bio.Phylo import TreeConstruction
except ImportError:
    HAVE_BIOPYTHON = False

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
        kwargs:
            can have: fusion
            most have: either conn_string
                       or all of these: "host", "dbname", "user", "password"

        """

        self.conn = None
        self.cur = None
        self.connstring = None
        self.fusion_url = None

        if kwargs.has_key('fusion') == True:
            self.fusion_url = kwargs['fusion']

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

        try:
            self.conn = psycopg2.connect(self.connstring)
            self.cur = self.conn.cursor(cursor_factory=DictCursor)
        except psycopg2.OperationalError as ex:
            raise SnapperDBInterrogationError(str(ex))

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
        self.cur.execute(sql, (sam_name, ))
        if self.cur.rowcount < 1:
            raise SnapperDBInterrogationError("No clustering information found for sample %s"
                                              % (sam_name))
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
            if self.fusion_url == None:
                distances = get_relevant_distances(self.cur, samid)
            else:
                distances = get_relevant_distances_fusion(self.cur, samid, self.fusion_url)
            sql = "SELECT pk_id, sample_name FROM samples"
            self.cur.execute(sql)
            id2name = {r['pk_id']: r['sample_name'] for r in self.cur.fetchall()}
        else:
            if self.fusion_url == None:
                distances = get_distances(self.cur, samid, list(close_samples))
            else:
                distances = get_distances_fusion(samid, list(close_samples), self.fusion_url)
        result_samples = distances[:neighbours]

        for (sa, di) in distances[neighbours:]:
            if di == result_samples[-1][1]:
                result_samples.append((sa, di))

        result_samples = [(id2name[sa], di) for (sa, di) in result_samples if sa != samid]
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
        self.cur.execute(sql, (sam_name, ))
        if self.cur.rowcount < 1:
            raise SnapperDBInterrogationError("No clustering information found for sample %s"
                                              % (sam_name))
        row = self.cur.fetchone()
        snad = [row['t0'], row['t5'], row['t10'], row['t25'], row['t50'], row['t100'], row['t250']]
        samid = row['pk_id']

        ct = get_closest_threshold(dis)
        if ct != None:
            # selected distance <250 -> use only samples in associated cluster for calculation
            t_ct = 't%i' % (ct)
            cluster = snad[levels.index(ct)]
            sql = "SELECT s.sample_name AS samname, c.fk_sample_id AS samid FROM sample_clusters c, samples s WHERE c."+t_ct+"=%s AND s.pk_id=c.fk_sample_id"
            self.cur.execute(sql, (cluster, ))
        else:
            # selected distance >250 -> use all samples that have been clustered and are not ignored for calculation
            sql = "SELECT s.sample_name AS samname, c.fk_sample_id AS samid FROM sample_clusters c, samples s WHERE s.pk_id=c.fk_sample_id AND s.ignore_sample IS FALSE AND s.sample_name<>%s"
            self.cur.execute(sql, (sam_name, ))

        id2name = {}
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
            logging.info("Calculating distances to %i samples.", len(neighbours))
            if self.fusion_url == None:
                distances = get_distances(self.cur, samid, neighbours)
            else:
                distances = get_distances_fusion(samid, neighbours, self.fusion_url)
            result_samples = [(id2name[s], d) for (s, d) in distances if d <= dis]
            if len(result_samples) <= 0:
                logging.info("No samples found this close to the query sample.")
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
            "1.2.3.4.5.6.7" if successful else None
        """

        sql = "SELECT c.t0, c.t5, c.t10, c.t25, c.t50, c.t100, c.t250 FROM sample_clusters c, samples s WHERE s.pk_id=c.fk_sample_id AND s.sample_name=%s"
        self.cur.execute(sql, (sam_name, ))
        if self.cur.rowcount < 1:
            raise SnapperDBInterrogationError("No clustering information found for sample %s"
                                              % (sam_name))
        elif self.cur.rowcount > 1:
            raise SnapperDBInterrogationError("Too much clustering information found for sample %s"
                                              % (sam_name))
        else:
            row = self.cur.fetchone()
            return "%i.%i.%i.%i.%i.%i.%i" % (row['t250'],
                                             row['t100'],
                                             row['t50'],
                                             row['t25'],
                                             row['t10'],
                                             row['t5'],
                                             row['t0'])

    # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

    def get_nearest(self, sam_name):
        """
        Derive the nearest distance for this sample from the SNP address for the sample.

        Parameters
        ----------
        sam_name: str
            name of the query sample

        Returns
        -------
        nearest: str
            e.g.: 100>=x>50
        """

        sql = "SELECT c.t0, c.t5, c.t10, c.t25, c.t50, c.t100, c.t250 FROM sample_clusters c, samples s WHERE s.pk_id=c.fk_sample_id AND s.sample_name=%s"
        self.cur.execute(sql, (sam_name, ))
        if self.cur.rowcount < 1:
            raise SnapperDBInterrogationError("No clustering information found for sample %s"
                                              % (sam_name))
        elif self.cur.rowcount > 1:
            raise SnapperDBInterrogationError("Too much clustering information found for sample %s"
                                              % (sam_name))
        else:
            row = self.cur.fetchone()

        levels = [250, 100, 50, 25, 10, 5, 0]
        for lvl in levels:
            t_lvl = 't%i' % (lvl)
            sql = "SELECT pk_id FROM sample_clusters WHERE " + t_lvl + "=%s"
            self.cur.execute(sql, (row[t_lvl], ))
            if self.cur.rowcount == 1:
                if lvl == 250:
                    nearest = "x>250"
                else:
                    nearest = "%i>=x>%i" % (levels[levels.index(lvl)-1], lvl)
                break
        # yes, for-loop-else, will be executed if the for loop did not exit via the *break*
        else:
            nearest = "x=0"

        return nearest

    # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

    def get_sample_history(self, sam_name):
        """
        Get the current snp address for a sample and any possible previous snp addresses.

        Parameters
        ----------
        sam_name: str
            name of the query sample

        Returns
        -------
        res: dict
            {'current_snad': 1.2.3.4.5.6.7,
             'history': [{'old': '100.2.3.4.5.6.99', 'new': '1.2.3.4.5.6.99',
                          'time': 2017-09-22 15:56:22.427083},
                         {'old': '1.2.3.4.5.6.99',   'new': '1.2.3.4.5.6.7',
                          'time': 2017-09-23 12:00:22.427083},
                         ...]}
        """

        sql = "SELECT s.pk_id, c.t0, c.t5, c.t10, c.t25, c.t50, c.t100, c.t250 \
               FROM sample_clusters c, samples s WHERE s.pk_id=c.fk_sample_id AND s.sample_name=%s"
        self.cur.execute(sql, (sam_name, ))
        if self.cur.rowcount < 1:
            raise SnapperDBInterrogationError("No clustering information found for sample %s"
                                              % (sam_name))
        elif self.cur.rowcount > 1:
            raise SnapperDBInterrogationError("Too much clustering information found for sample %s"
                                              % (sam_name))
        else:
            row = self.cur.fetchone()

        sam_id = row['pk_id']
        levels = [250, 100, 50, 25, 10, 5, 0]
        res = {'current_snad': '.'.join([str(row['t%i' % (lvl)]) for lvl in levels]),
               'history': []}

        sql = "SELECT t0_old, t5_old, t10_old, t25_old, t50_old, t100_old, t250_old, t0_new, \
                      t5_new, t10_new, t25_new, t50_new, t100_new, t250_new, renamed_at \
                      FROM sample_history WHERE fk_sample_id=%s"
        self.cur.execute(sql, (sam_id, ))
        rows = self.cur.fetchall()
        for r in rows:
            res['history'].append({'old': '.'.join([str(r['t%i_old' % (lvl)]) for lvl in levels]),
                                   'new': '.'.join([str(r['t%i_new' % (lvl)]) for lvl in levels]),
                                   'time': r['renamed_at']})

        return res

    # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

    def get_tree(self, samples, clusters, method, **kwargs):
        """
        Build a tree for a given set of samples with a given method.

        Parameters
        ----------
        samples: list
            list of sample names
        clusters: dict
            {'lvl': [1,2,3], 'lvl': [2,5,8]}
        method: str
            either 'ML' or 'NJ'

        Returns
        -------
        tree as newick string

        """

        treesams = {}

        if samples != None:

            sql = "SELECT pk_id, sample_name FROM samples WHERE sample_name IN %s"
            self.cur.execute(sql, (tuple(samples), ))
            rows = self.cur.fetchall()
            for r in rows:
                treesams[r['sample_name']] = r['pk_id']

            missing = set(samples).difference(set(treesams.keys()))
            if len(missing) > 0:
                logging.warning("The following sample names were not found in the database: %s",
                                str(list(missing)))
            else:
                logging.info("All samples names provided were found in the database.")

        if clusters != None:
            for t_lvl, clusterlist in clusters.items():
                try:
                    sql = "SELECT c.fk_sample_id AS id, s.sample_name AS name FROM sample_clusters c, samples s WHERE s.pk_id=c.fk_sample_id AND c."+t_lvl+" IN %s"
                    self.cur.execute(sql, (tuple(clusterlist), ))
                    rows = self.cur.fetchall()
                    for r in rows:
                        treesams[r['name']] = r['id']
                except psycopg2.ProgrammingError as e:
                    raise SnapperDBInterrogationError(e)

        nofsams = len(treesams.keys())
        if nofsams < 3:
            raise SnapperDBInterrogationError("At least 3 samples are required to make a tree. Only %i found." % nofsams)
        elif nofsams > 400:
            if kwargs.has_key('overwrite_max') == True and kwargs['overwrite_max'] == True:
                pass
            else:
                raise SnapperDBInterrogationError("This tree would contain %i samples. A maximum of 400 is permitted. Please select a more targeted subset." % nofsams)
        else:
            pass

        logging.info("Calculating %s tree for %i samples.", method, nofsams)

        if method == 'NJ':
            if HAVE_BIOPYTHON == False:
                raise SnapperDBInterrogationError("You need to have Biopython for making NJ trees.")
            return self._make_nj_tree(treesams, kwargs['dm'])
        elif method == 'ML':
            if self._can_we_make_an_ml_tree() == False:
                raise SnapperDBInterrogationError("You need to have FastTree for making ML trees.")
            return self._make_ml_tree(treesams, kwargs['ref'], kwargs['refname'], kwargs['rmref'])
        else:
            raise SnapperDBInterrogationError("%s is an unsupported method." % (method))

    # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

    def _make_nj_tree(self, treesams, dm):
        """
        **PRIVATE**

        Parameters
        ----------
        treesams: dict
            {sam name: samid, sam name: samid, ...}

        Returns
        -------
        nwkstring: str
            tree as newick string
        """

        iNofSams = len(treesams.keys())
        logging.info("Calculating %i distances. Patience!", ((iNofSams**2) - iNofSams) / 2)

        dist_mat = None
        if self.fusion_url == None:
            dist_mat = get_distance_matrix(self.cur, treesams.values())
        else:
            dist_mat = get_distance_matrix_fusion(treesams.values(), self.fusion_url)

        if dm != None:
            logging.info("Distance matrix written to file: %s", dm)
            if os.path.exists(dm) == True:
                os.remove(dm)

        aSampleNames = treesams.keys()
        aSimpleMatrix = []
        for i, sample_1 in enumerate(aSampleNames):
            mat_line = []
            for j, sample_2 in enumerate(aSampleNames):
                if j < i:
                    sid1 = treesams[sample_1]
                    sid2 = treesams[sample_2]
                    mat_line.append(dist_mat[sid1][sid2])
                elif j == i:
                    mat_line.append(0)
                else:
                    pass
            aSimpleMatrix.append(mat_line)
            if dm != None:
                with open(dm, 'a') as f:
                    f.write("%s\n" % ','.join([sample_1] + [str(x) for x in mat_line[:-1]]))

        logging.info("Bulding tree.")
        oDistMat = TreeConstruction._DistanceMatrix(aSampleNames, aSimpleMatrix)
        constructor = TreeConstruction.DistanceTreeConstructor()
        oTree = constructor.nj(oDistMat)

        # I don't know how to get newick string from this object without a file ...
        td = tempfile.mkdtemp()
        tmpfile = os.path.join(td, 'tree.nwk')
        Phylo.write(oTree, tmpfile, 'newick')
        nwkstring = ""
        with open(tmpfile, 'r') as f:
            nwkstring = f.read()
        shutil.rmtree(td)

        return nwkstring

    # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

    def _make_ml_tree(self, treesams, ref, refname, rmref):
        """
        **PRIVATE**

        Parameters
        ----------
        treesams: dict
            {samid: sam name, samid: sam name, ...}

        Returns
        -------

        not sure yet

        """

        td = tempfile.mkdtemp()
        tmpfile = os.path.join(td, 'align.fasta')

        inp = {'whole_genome': False,
               'snp_address': False,
               'remove_invariant_npos': True,
               'reference': ref,
               'column_gaps': 0.5,
               'column_Ns': 0.5,
               'cmd': 'get_alignment',
               'db': self.connstring,
               'sample_Ns': None,
               'sample_gaps': None,
               'version': 'internal',
               'include': None,
               'samples': treesams.keys(),
               'exclude': None,
               'debug': False,
               'sample_Ns_gaps_auto_factor': 2.0,
               'name_of_ref_in_db': refname,
               'out': tmpfile,
               'remove_ref': rmref}

        if get_alignment.main(inp) != 0:
            raise SnapperDBInterrogationError("Error in get_alignment main.")

        cmd = "FastTree -nt %s" % (tmpfile)
        logging.info("Running FastTree now. Patience.")
        p = subprocess.Popen(cmd, shell=True, stdin=None,
                             stdout=subprocess.PIPE,
                             stderr=subprocess.PIPE, close_fds=True)
        (p_out, p_err) = p.communicate()

        logging.debug("FastTree stderr: %s", p_err)

        shutil.rmtree(td)

        return p_out

    # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

    def _can_we_make_an_ml_tree(self):
        """
        **PRIVATE**

        Check if we can call FastTree.

        Parameters
        ----------

        Returns
        -------
        true if we can, else false

        """

        p = subprocess.Popen("which FastTree", shell=True, stdin=None,
                             stdout=subprocess.PIPE,
                             stderr=subprocess.PIPE, close_fds=True)
        (p_out, p_err) = p.communicate()

        flag = not "no FastTree in" in p_err

        if flag == True:
            logging.info("Using %s to make an ML tree.", p_out.strip())

        return flag

# --------------------------------------------------------------------------------------------------
