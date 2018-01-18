-- table definitions for all snapperdb3 databases

DROP TABLE IF EXISTS samples;
CREATE TABLE samples (
    pk_id SERIAL PRIMARY KEY,
    sample_name text,
    molis_id text,
    ngs_id integer,
    ignore_sample boolean DEFAULT FALSE,
    ignore_zscore boolean DEFAULT FALSE,
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
CREATE TABLE sample_clusters (
    pk_id SERIAL PRIMARY KEY,
    fk_sample_id integer references samples(pk_id),
    t0 integer,
    t5 integer,
    t10 integer,
    t25 integer,
    t50 integer,
    t100 integer,
    t250 integer,
    t0_mean double precision,
    t5_mean double precision,
    t10_mean double precision,
    t25_mean double precision,
    t50_mean double precision,
    t100_mean double precision,
    t250_mean double precision
);

DROP TABLE IF EXISTS cluster_stats;
CREATE TABLE cluster_stats (
    pk_id SERIAL PRIMARY KEY,
    cluster_level text,
    cluster_name integer,
    nof_members integer,
    nof_pairwise_dists integer,
    mean_pwise_dist double precision,
    stddev double precision
);

DROP TABLE IF EXISTS merge_log;
CREATE TABLE merge_log (
    pk_id SERIAL PRIMARY KEY,
    cluster_level text,
    source_cluster integer,
    target_cluster integer,
    time_of_merge timestamp
);

DROP TABLE IF EXISTS sample_history;
CREATE TABLE sample_history (
    pk_id SERIAL PRIMARY KEY,
    fk_sample_id integer references samples(pk_id),
    t0_old integer,
    t5_old integer,
    t10_old integer,
    t25_old integer,
    t50_old integer,
    t100_old integer,
    t250_old integer,
    t0_new integer,
    t5_new integer,
    t10_new integer,
    t25_new integer,
    t50_new integer,
    t100_new integer,
    t250_new integer,
    renamed_at timestamp
);

-- #################################################################################################
-- adding trees table

DROP TABLE IF EXISTS trees;
CREATE TABLE trees (
    pk_id SERIAL PRIMARY KEY,
    nwkfile bytea,
	t5_name integer,
	sample_set integer[],
	mod_date timestamp,
	created_at timestamp,
	lockdown boolean DEFAULT FALSE
);


CREATE EXTENSION intarray;
