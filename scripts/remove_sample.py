import logging

import argparse
from argparse import RawTextHelpFormatter
import psycopg2
from psycopg2.extras import DictCursor

from lib.distances import get_distances
from lib.ClusterStats import ClusterStats

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

    return r'''This removes a samples from the database, which is not a trivial thing to do because
cluster stats and stats for other samples need to be updated.'''

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

    args = argparse.ArgumentParser(description=get_desc(), formatter_class=RawTextHelpFormatter)

    args.add_argument("--connstring",
                      "-c",
                      type=str,
                      metavar="CONNECTION",
                      required=True,
                      dest="db",
                      help="Connection string for db. REQUIRED")

    args.add_argument("--sample",
                      "-s",
                      type=str,
                      required=True,
                      dest="sample",
                      help="Name of sample to remove. REQUIRED.")

    args.add_argument("--just-ignore",
                      action='store_true',
                      help="""Sample and variant information will be retained in database,
but clustering information will be removed. Ignore_sample will
be set to TRUE. [DEFAULT: Remove everything. Sample can be added
and clustered again later.]""")

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

    try:
        # open db
        conn = psycopg2.connect(args['db'])
        cur = conn.cursor(cursor_factory=psycopg2.extras.DictCursor)

        # get the pkid of this sample
        sql = "SELECT pk_id, ignore_sample, ignore_zscore FROM samples WHERE sample_name=%s"
        cur.execute(sql, (args['sample'], ))
        if cur.rowcount < 1:
            logging.error("Sample '%s' not found in database.", args['sample'])
            return 1
        row = cur.fetchone()

        sample_id = row['pk_id']
        igno_flag = row['ignore_sample']
        zscr_flag = row['ignore_zscore']

        # is this sample is already ignored
        if igno_flag == True:
            logging.info("This sample is already ignored.")
            if args['just_ignore'] == True:
                logging.info("You chose not to remove it completely.")
                return 0
            else:
                drop_sample(cur, sample_id)
                conn.commit()
                return 0


        if zscr_flag == False:
            # get the snp address
            sql = "SELECT t0, t5, t10, t25, t50, t100, t250 FROM sample_clusters WHERE fk_sample_id=%s"
            cur.execute(sql, (sample_id, ))
            if cur.rowcount != 1:
                logging.error("There is not exactly one entry in sample clusters for this samples. :-(")
                return 1
            row = cur.fetchone()
            snad = [row['t0'], row['t5'], row['t10'], row['t25'], row['t50'], row['t100'], row['t250']]


            sql = "SELECT c.fk_sample_id AS samid FROM sample_clusters c, samples s WHERE c.t250=%s AND c.fk_sample_id=s.pk_id AND s.ignore_zscore IS FALSE"
            cur.execute(sql, (row['t250'], ))
            t250_members = [r['samid'] for r in cur.fetchall()]

            try:
                t250_members.remove(sample_id)
            except ValueError:
                logging.error("Bizzare data inconsistency for sample id %s and t250 cluster %s.", sample_id, row['t250'])
                return 1

            logging.info("Calculating %i distances to update stats.", len(t250_members))
            distances = get_distances(cur, sample_id, t250_members)

            # update all stats in all clusters on all levels
            for clu, lvl in zip(snad, [0, 5, 10, 25, 50, 100, 250]):
                if update_cluster_stats_post_removal(cur, sample_id, clu, lvl, distances) == None:
                    logging.error("Problem with updating cluster stats.")
                    return None
        else:
            logging.info("This sample previously failed the zscore check and is not considered in stats. Skipping stats update.")

        # now remove the sample
        sql = "DELETE FROM sample_clusters WHERE fk_sample_id=%s"
        cur.execute(sql, (sample_id, ))

        if args['just_ignore'] == True:
            logging.info("You chose not to remove the sample completely.")
            sql = "UPDATE samples SET ignore_sample=True WHERE pk_id=%s"
            cur.execute(sql, (sample_id, ))
        else:
            drop_sample(cur, sample_id)

        conn.commit()

    except psycopg2.Error as e:
         logging.error("Database reported error: %s" % (str(e)))
         return None
    finally:
        # close all dbs
        cur.close()
        conn.close()

    return 0

# end of main --------------------------------------------------------------------------------------

def update_cluster_stats_post_removal(cur, sid, clu, lvl, distances):
    """
    Update the cluster stats and the sample stats for removing the sample from the cluster.

    Parameters
    ----------
    cur: obj
        database cursor
    sid: int
        pk id of sample to remove
    clu: int
        cluster id to remove from
    lvl: int
        cluster level of removal
    distances: list of tuples
        [(samid, distance), (samid, distance), (samid, distance), ...]

    Returns
    -------
    0 if fine
    None if fail
    """

    t_lvl = "t%i" % (lvl)

    logging.info("Updating stats for cluster %s on level %s.", clu, t_lvl)

    # get the cluster stats from the database
    sql = "SELECT nof_members, mean_pwise_dist, stddev FROM cluster_stats WHERE cluster_level=%s AND cluster_name=%s"
    cur.execute(sql, (t_lvl, clu, ))
    if cur.rowcount == 0:
        logging.error("Cluster stats for level %s and cluster %s not found.", t_lvl, clu)
        return None
    row = cur.fetchone()

    # when deleting the last member of this cluster there is no need to update anything, just get rid of it
    if row['nof_members'] <= 1:
        logging.debug("This is the last member of cluster %s on level %s. Deleting cluster stats.", clu, t_lvl)
        sql = "DELETE FROM cluster_stats WHERE cluster_level=%s AND cluster_name=%s"
        cur.execute(sql, (t_lvl, clu, ))
        return 0

    # create cluster stats object from the information in the database
    oStats = ClusterStats(members=row['nof_members'], stddev=row['stddev'], mean=row['mean_pwise_dist'])

    # get other members of this cluster, we know it must be at least one
    sql = "SELECT c.fk_sample_id AS samid FROM sample_clusters c, samples s WHERE c."+t_lvl+"=%s AND c.fk_sample_id=s.pk_id AND s.ignore_zscore IS FALSE"
    cur.execute(sql, (clu, ))
    members = [r['samid'] for r in cur.fetchall()]
    try:
        members.remove(sid)
    except ValueError:
        logging.error("Bizzare data inconsistency for sample id %s and %s cluster %s.", sid, t_lvl, clu)
        return None

    # get all distances fromt he sample to be removed to all other members of the cluster
    # ad update the stats object with this information
    this_di = [d for (s, d) in distances if s in members]
    assert oStats.members == (len(this_di) + 1)
    oStats.remove_member(this_di)

    # then update the cluster stats in the database with the info from the object
    sql = "UPDATE cluster_stats SET (nof_members, nof_pairwise_dists, mean_pwise_dist, stddev) = (%s, %s, %s, %s) WHERE cluster_level=%s AND cluster_name=%s"
    cur.execute(sql, (oStats.members, oStats.nof_pw_dists, oStats.mean_pw_dist, oStats.stddev_pw_dist, t_lvl, clu, ))

    # for all othr members of this cluster
    for mem in members:

        # get the mean distance to all other members
        sql = "SELECT "+t_lvl+"_mean FROM sample_clusters WHERE fk_sample_id=%s"
        cur.execute(sql, (mem, ))
        if cur.rowcount == 0:
            logging.error("Cluster %s not found in sample_clusters table.", mem)
            return None
        row = cur.fetchone()
        p_mean = row[t_lvl+'_mean']

        # update this mean by removing one distance and update the database table
        x = [d for (s, d) in distances if s == mem][0]
        try:
            n_mean = ((p_mean * len(members)) - x) / float(len(members) - 1)
        except ZeroDivisionError:
            n_mean = None

        sql = "UPDATE sample_clusters SET "+t_lvl+"_mean=%s WHERE fk_sample_id=%s"
        cur.execute(sql, (n_mean, mem, ))

    return 0

# --------------------------------------------------------------------------------------------------

def drop_sample(cur, sid):
    """
    Remove the sample with the id from the variants and samples tables.

    Parameters
    ----------
    cur: obj
        database cursor
    sid: int
        pk id of sample to remove

    Returns
    -------
    0
    """

    logging.info("Removing the sample from the samples and variants table.")

    sql = "DELETE FROM variants WHERE fk_sample_id=%s"
    cur.execute(sql, (sid, ))

    sql = "DELETE FROM samples WHERE pk_id=%s"
    cur.execute(sql, (sid, ))

    return 0

# --------------------------------------------------------------------------------------------------

if __name__ == "__main__":
    exit(main(vars(get_args().parse_args())))
