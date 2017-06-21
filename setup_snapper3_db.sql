-- table definitions for all snapperdb3 databases

DROP TABLE IF EXISTS samples;
CREATE TABLE samples (
    pk_id SERIAL PRIMARY KEY,
    sample_name text,
    molis_id text,
    ngs_id integer,
    ignore_sample boolean DEFAULT FALSE,
    date_added timestamp
);

DROP TABLE IF EXISTS contigs;
CREATE TABLE contigs (
    pk_id SERIAL PRIMARY KEY,
    name text,
    length integer
);

DROP TABLE IF EXISTS variants;
CREATE TABLE variants (
    pk_id SERIAL PRIMARY KEY,
    fk_sample_id integer references samples(pk_id),
    fk_contig_id integer references contigs(pk_id),
    a_pos integer[],
    c_pos integer[],
    g_pos integer[],
    t_pos integer[],
    n_pos integer[],
    gap_pos integer[]
);

DROP TABLE IF EXISTS sample_clusters;
CREATE TABLE sample_clusters  (
    pk_id SERIAL PRIMARY KEY,
    fk_sample_id integer references samples(pk_id),
    t0 integer,
    t5 integer,
    t10 integer,
    t25 integer,
    t50 integer,
    t100 integer,
    t250 integer,
    t0_mean real,
    t5_mean real,
    t10_mean real,
    t25_mean real,
    t50_mean real,
    t100_mean real,
    t250_mean real
);

DROP TABLE IF EXISTS cluster_stats;
CREATE TABLE cluster_stats (
    pk_id SERIAL PRIMARY KEY,
    cluster_level text,
    cluster_name integer,
    nof_members integer,
    nof_pairwise_dists integer,
    mean_pwise_dist real,
    stddev real
);

CREATE EXTENSION intarray;
