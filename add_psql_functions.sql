-- #################################################################################################
-- array_symdiff helper

DROP FUNCTION IF EXISTS array_symdiff(anyarray, anyarray);
CREATE OR REPLACE FUNCTION public.array_symdiff(
    anyarray,
    anyarray)
  RETURNS anyarray AS
$BODY$

    SELECT (
        ($1 | $2)
        -
        ($1 & $2)
    );

$BODY$
  LANGUAGE sql
  COST 100;

-- #################################################################################################
-- get distance from one sample to a list of others by submitting variants in list format

DROP FUNCTION IF EXISTS public.get_sample_distances_by_lists(
    IN in_a integer[],
    IN in_c integer[],
    IN in_g integer[],
    IN in_t integer[],
    IN in_n integer[],
    IN in_gap integer[],
    IN in_chr_id integer,
    IN in_samples integer[],
    OUT integer,
    OUT integer,
    OUT integer);
CREATE OR REPLACE FUNCTION public.get_sample_distances_by_lists(
    IN in_a integer[],
    IN in_c integer[],
    IN in_g integer[],
    IN in_t integer[],
    IN in_n integer[],
    IN in_gap integer[],
    IN in_chr_id integer,
    IN in_samples integer[],
    OUT integer,
    OUT integer,
    OUT integer)
  RETURNS SETOF record AS
$BODY$
BEGIN
 RETURN QUERY WITH summary AS (SELECT fk_sample_id,
		fk_contig_id,
		a_pos,
		c_pos,
		g_pos,
		t_pos,
		n_pos,
		gap_pos,
		array_length(
		(
			array_symdiff(in_a, a_pos)
			|
			array_symdiff(in_c, c_pos)
			|
			array_symdiff(in_g, g_pos)
			|
			array_symdiff(in_t, t_pos)
		)
		-
		(in_n | n_pos | gap_pos | in_gap)
	,1) AS dist FROM variants WHERE variants.fk_contig_id=in_chr_id AND variants.fk_sample_id= ANY (in_samples))
SELECT
	fk_sample_id,
	fk_contig_id,
	dist
FROM
	summary
ORDER BY
	dist ASC;
END;
$BODY$
  LANGUAGE plpgsql IMMUTABLE
  COST 100
  ROWS 1000;

-- #################################################################################################
-- get distance from one sample to a list of others by submitting the sample id

DROP FUNCTION IF EXISTS public.get_sample_distances_by_id(
    IN pivot integer,
    IN in_chr_id integer,
    IN in_samples integer[],
    OUT integer,
    OUT integer,
    OUT integer);
CREATE OR REPLACE FUNCTION public.get_sample_distances_by_id(
    IN pivot integer,
    IN in_chr_id integer,
    IN in_samples integer[],
    OUT integer,
    OUT integer,
    OUT integer)
  RETURNS SETOF record AS
$BODY$
BEGIN
 RETURN QUERY WITH
getone AS (SELECT
            a_pos AS in_a,
	        c_pos AS in_c,
	        g_pos AS in_g,
	        t_pos AS in_t,
	        n_pos AS in_n,
	        gap_pos AS in_gap
	        FROM variants WHERE fk_sample_id=pivot AND fk_contig_id=in_chr_id)
SELECT
        fk_sample_id,
		fk_contig_id,
		array_length(
		(
			array_symdiff(in_a, a_pos)
			|
			array_symdiff(in_c, c_pos)
			|
			array_symdiff(in_g, g_pos)
			|
			array_symdiff(in_t, t_pos)
		)
		-
		(in_n | n_pos | gap_pos | in_gap)
	,1) AS dist FROM variants, getone
	WHERE variants.fk_contig_id=in_chr_id
	    AND variants.fk_sample_id=ANY(in_samples)
	    ORDER BY dist ASC;
END;
$BODY$
  LANGUAGE plpgsql IMMUTABLE
  COST 100
  ROWS 1000;

-- #################################################################################################
-- get pairwise distance for two sample ids

DROP FUNCTION IF EXISTS public.get_pairwise_distance(
    IN in_chr_id integer,
    IN in_sample1 integer,
    IN in_sample2 integer,
    OUT integer);
CREATE OR REPLACE FUNCTION public.get_pairwise_distance(
    IN in_chr_id integer,
    IN in_sample1 integer,
    IN in_sample2 integer,
    OUT integer)
  RETURNS integer AS
$BODY$
 WITH samqryone AS
	(SELECT a_pos AS apo,
	        c_pos AS cpo,
	        g_pos AS gpo,
	        t_pos AS tpo,
	        n_pos AS npo,
	        gap_pos AS gappo
	        FROM variants WHERE fk_sample_id=$2 AND fk_contig_id=$1),
	samqrytwo AS
	(SELECT a_pos AS apt,
	        c_pos AS cpt,
	        g_pos AS gpt,
	        t_pos AS tpt,
	        n_pos AS npt,
	        gap_pos AS gappt
	        FROM variants WHERE fk_sample_id=$3 AND fk_contig_id=$1)
SELECT array_length(
		(array_symdiff(apo, apt)
		| array_symdiff(cpo, cpt)
		| array_symdiff(gpo, gpt)
		| array_symdiff(tpo, tpt))
		 - (npo | gappo | npt | gappt), 1)
		FROM samqryone, samqrytwo;
$BODY$
  LANGUAGE sql IMMUTABLE
  COST 100;

-- #################################################################################################
-- get distance from one sample to a list of others by submitting the sample id

DROP FUNCTION IF EXISTS public.get_all_distances_by_id(
    IN pivot integer,
    IN in_chr_id integer,
    OUT integer,
    OUT integer,
    OUT integer);
CREATE OR REPLACE FUNCTION public.get_all_distances_by_id(
    IN pivot integer,
    IN in_chr_id integer,
    OUT integer,
    OUT integer,
    OUT integer)
  RETURNS SETOF record AS
$BODY$
BEGIN
 RETURN QUERY WITH
getone AS (SELECT
            a_pos AS in_a,
	        c_pos AS in_c,
	        g_pos AS in_g,
	        t_pos AS in_t,
	        n_pos AS in_n,
	        gap_pos AS in_gap
	        FROM variants WHERE fk_sample_id=pivot AND fk_contig_id=in_chr_id)
SELECT
        fk_sample_id,
		fk_contig_id,
		array_length(
		(
			array_symdiff(in_a, a_pos)
			|
			array_symdiff(in_c, c_pos)
			|
			array_symdiff(in_g, g_pos)
			|
			array_symdiff(in_t, t_pos)
		)
		-
		(in_n | n_pos | gap_pos | in_gap)
	,1) AS dist FROM variants, getone
	WHERE variants.fk_contig_id=in_chr_id and fk_sample_id!=pivot ORDER BY dist ASC;
END;
$BODY$
  LANGUAGE plpgsql IMMUTABLE
  COST 100
  ROWS 1000;

-- #################################################################################################
