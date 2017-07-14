"""
File contains some to get some distance calculations done for snapperdb v3.

author: ulf.schaefer@phe.gov.uk

"""

from lib.ClusterStats import ClusterStats

# --------------------------------------------------------------------------------------------------

def register_sample(cur, sample_id, distances, new_snad, levels=[0, 5, 10, 25, 50, 100, 250]):
    """
    Registers a sample in the database and updates all cluster and sample statistics

    Parameters
    ----------
    cur: obj
        database cursor
    sample_id: int
        pk_id in samples tables
    distances: list of tuples
        sorted list of tuples with (sample_id, distance) with closes sample first
        e.g. [(298, 0), (37, 3), (55, 4)]
    new_snad: list of 7 ints
        preliminary snp address

    Returns
    -------
    final_snad: list of 7 ints
        the final snp address of the sample
    None if problem
    """

    final_snad = {}
    means = {}

    for lvl, cluster in zip(levels, new_snad):
        print lvl, cluster

        t_lvl = 't%i' % (lvl)

        if cluster == None:
            sql = "SELECT max("+t_lvl+") AS m FROM sample_clusters"
            cur.execute(sql)
            row = cur.fetchone()
            final_snad[lvl] = row['m'] + 1
            means[lvl] = None
        else:
            # get current cluster stats
            sql = "SELECT nof_members, mean_pwise_dist, stddev FROM cluster_stats WHERE cluster_name=%s AND cluster_level=%s"
            cur.execute(sql, (cluster, t_lvl))
            if cur.rowcount != 1:
                logging.error("Uncertain about stats for %s level %s", cluster, t_lvl)
                return None
            statsrow = cur.fetchone()

            nof_members = statsrow['nof_members']

            # get the current members of the cluster and the distances to them
            sql = "SELECT c.fk_sample_id FROM sample_clusters c, samples s WHERE c."+t_lvl+"=%s AND s.pk_id=c.fk_sample_id AND s.ignore_zscore IS FALSE"
            cur.execute(sql, (cluster, ))
            rows = cur.fetchall()
            current_members = [r['fk_sample_id'] for r in rows]
            dis_to_cu_mems = [d for (s, d) in distances if s in current_members]

            if nof_members > 1:
                # create cluster stats object
                oStats = ClusterStats(members=nof_members, stddev=statsrow['stddev'], mean=statsrow['mean_pwise_dist'])
            else:
                oStats = ClusterStats(members=nof_members, dists=[])

            # update the cluster stats object
            oStats.add_member(dis_to_cu_mems)

            # update cluster stats table
            sql  = "UPDATE cluster_stats SET (nof_members, nof_pairwise_dists, mean_pwise_dist, stddev) = (%s, %s, %s, %s) WHERE cluster_name=%s AND cluster_level=%s"
            cur.execute(sql, (oStats.members, oStats.nof_pw_dists, oStats.mean_pw_dist, oStats.stddev_pw_dist, cluster, t_lvl, ))

            means[lvl] = sum(dis_to_cu_mems) / float(len(dis_to_cu_mems))
            final_snad[lvl] = cluster

            # update the stats of all other members of the cluster
            for o_mem in current_members:

                sql = "SELECT "+t_lvl+"_mean FROM sample_clusters WHERE fk_sample_id=%s"
                cur.execute(sql, (o_mem, ))
                if cur.rowcount != 1:
                    logging.error("Uncertain about clustering info for sample %s", cluster, t_lvl)
                    return None
                row = cur.fetchone()
                old_mean = row[t_lvl+'_mean']

                # old mean is None when the cluster has only one member
                if old_mean == None:
                    old_mean = 0.0

                new_dis = [d for (s, d) in distances if s == o_mem][0]
                new_mean = ((old_mean * (nof_members - 1)) + new_dis) / float(nof_members)

                sql = "UPDATE sample_clusters SET "+t_lvl+"_mean=%s WHERE fk_sample_id=%s"
                cur.execute(sql, (new_mean, o_mem, ))

    # end for lvl, cluster in zip(levels, new_snad):

    sql = "INSERT INTO sample_clusters (fk_sample_id, t0, t5, t10, t25, t50, t100, t250, t0_mean, t5_mean, t10_mean, t25_mean, t50_mean, t100_mean, t250_mean) VALUES (%s, %s, %s, %s, %s, %s, %s, %s, %s, %s, %s, %s, %s, %s, %s)"
    cur.execute(sql, (sample_id,
                      final_snad[levels[0]],
                      final_snad[levels[1]],
                      final_snad[levels[2]],
                      final_snad[levels[3]],
                      final_snad[levels[4]],
                      final_snad[levels[5]],
                      final_snad[levels[6]],
                      means[levels[0]],
                      means[levels[1]],
                      means[levels[2]],
                      means[levels[3]],
                      means[levels[4]],
                      means[levels[5]],
                      means[levels[6]],))

    return [final_snad[x] for x in levels]

# --------------------------------------------------------------------------------------------------
