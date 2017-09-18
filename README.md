# snapper3

A beta version of the software is provided in **/phengs/hpc_software/phe/snapper3/3-0**. To use it:

```bash
module load phe/snapper3/3-0
```

Prerequisites/dependencies:
- Python >= 2.7.6
- psycopg2 >= 2.5.2

The postgres server must be running **postgres >=9.6** and the contrib package must be installed so that the **intarray extension** can be created for each database.

All commands require a connection string of the format:

```bash
"host='db host IP' dbname='database_name' user='uname' password='password'"
```

New databases can be created with the script 'reset.sh'. The script can also be used to scrap existing databases and start from scratch.

```bash
sh reset.sh dbuser 158.119.123.123 my_snapper3_db
```

A script is provided to migrate an existing snapper v2 database to the new format.

```bash
usage: migrate_to_snapperV3.py [-h] --reference REFNAME --oldconnstring
                               CONNECTION --newconnstring CONNECTION

version 0.1, date 30Sep2016, author ulf.schaefer@phe.gov.uk

optional arguments:
  -h, --help            show this help message and exit
  --reference REFNAME, -r REFNAME
                        The sample_name of the reference genome in the
                        database.
  --oldconnstring CONNECTION, -o CONNECTION
                        Connection string for old db ('source')
  --newconnstring CONNECTION, -n CONNECTION
                        Connection string for new db ('target')

```

To create a new database, use the reset.sh script and then use the add_reference subcommand of snapper3.py

```bash
usage: snapper3.py add_reference [-h] --connstring CONNECTION --reference
                                 FASTAFILE --input JSONFILE [--ref-name NAME]

Takes variants for a sample in json format and adds them to the database.

optional arguments:
  -h, --help            show this help message and exit
  --connstring CONNECTION, -c CONNECTION
                        REQUIRED. Connection string for db.
  --reference FASTAFILE
                        REQUIRED. Fasta reference file.
  --input JSONFILE, -i JSONFILE
                        REQUIRED. Path to a input file.
  --ref-name NAME, -r NAME
                        The name of the reference to go into the db [default: reference file name before 1st dot]
```

The best way to add the variants of a new sample to the database is to submit them in json format. These json files can be made with the latest (1-4) version of Phenix from vcf files. The json file required here defines the ignore positions within the reference. It is usually made by calling phenix run_snp_pipeline with simulated reads on the reference and applying the usual filters to the vcf. This vcf is then converted to json and used here.

To add samples other than the reference to the database use:

```bash
usage: snapper3.py add_sample [-h] --input FILE --format FORMAT --connstring
                              CONNECTION --refname REFNAME
                              [--sample-name NAME] [--reference FASTAFILE]
                              [--min-coverage FLOAT]

Takes variants for a sample in json format and adds them to the database.

optional arguments:
  -h, --help            show this help message and exit
  --input FILE, -i FILE
                        REQUIRED. Path to a input file.
  --format FORMAT, -f FORMAT
                        REQUIRED. Choose from 'json' or 'fasta'.
  --connstring CONNECTION, -c CONNECTION
                        REQUIRED. Connection string for db.
  --refname REFNAME, -r REFNAME
                        REQUIRED. The sample_name of the reference genome in the database.
  --sample-name NAME, -s NAME
                        The name of the sample to go into the db [default: input file name before 1st dot]
  --reference FASTAFILE
                        Path to reference for this sample. Must be the same as used for the database.
                        REQUIRED when format is fasta, else ignored.
  --min-coverage FLOAT  Minimum coverage required to aloow sample in database. Only applicable with json
                        format, ignored for fasta.  This will check for the coverageMetaData annotation calculated by Phenix
                        and only allow the sample in, if the mean coverage is >= this value.
                        [default: do not check this]
```

Note:
- The name of the reference in the database is required here, because reference-ignore-positions are not stored for each sample individually.
- The sample can also be added from a fasta file (consensus sequence). This requires the reference as a fasta file (--reference).
- A minimum coverage can be specified (--min-coverage). This is read from the Phenix annotation in the json input file (-i). If the coverage is below the threshold the sample is not added to the database.

This will only create entries for the sample in the samples and variants table, but will do no clustering.

To cluster a sample use:

```bash
usage: snapper3.py cluster_sample [-h] --connstring CONNECTION --sample-name
                                  NAME [--no-zscore-check]
                                  [--with-registration] [--force-merge]

After the variants for a sample have been added to the database, use this to
determine the clustering for this sample. Will perform all statictical checks and
merging if necessary and update the database accordingly. If statistical checks fail
database will not be updated.

optional arguments:
  -h, --help            show this help message and exit
  --connstring CONNECTION, -c CONNECTION
                        Connection string for db. REQUIRED
  --sample-name NAME, -s NAME
                        The name of the sample to go into the db. REQUIRED.
  --no-zscore-check     Do not perform checks and just add the sample. It's fine.
                        [Default: Perform checks.]
  --with-registration   Register the clustering for this sample in the database
                        and update the cluster stats. [Default: Do not register.]
  --force-merge         Add the sample even if it causes clusters to merge.
                        [Default: Do not add if merge required.]
```

Note:
- **--no-zscore-check** will switch off all zscore calculations and cluster the sample. A sample added in this way will be excluded from all sample statistics and not be considered in future z-score checks.
- **--with-registration** without this option no permamant change will be made to the database. Effectively the postgres transaction is not commited before exiting.
- **--force-merge** The default behaviour is that if a sample requires a merge, the sample is not clustered. With this option the sample will be clustered and all merging will be done.

If a samples is dodgy and you want to not cluster it, use *remove_sample* to either remove all traces of the sample from the database or set the cluster to ignored=yes.

Further options:
- **get_alignment** Takes a list of sample names as a blanck separated list and provides an alignment in the same way that Phenix vcf2fasta does.
- **remove_sample** If you use this for a samples that has already been clustered, it removes the sample and un-does all the clustering. This **TAKES FOREVER** because it needs to be checked if removing this samples causes any clusters to be split. Do not use unless you have to (i.e. no backup available and staring from scratch not an option.)
- **export_sample_variants** If you want the variants for a sample in json format
- **get_closest** For a given sample, return the N closest samples in the database, or all samples cluster than x SNPs. This uses the SnapperDBInterrogation class.


API:
The same class also serves as the most glorious API. Use it like this:

```python
from lib.SnapperDBInterrogation import SnapperDBInterrogation, SnapperDBInterrogationError
with SnapperDBInterrogation(host=conf['db_host'],
                            dbname=dbname,
                            user=conf['db_username'],
                            password=conf['db_password']) as sdbi:
    print sdbi.get_closest_samples("123456_H1234567890-2", 10)
    print sdbi.get_samples_below_threshold("123456_H1234567890-2", 20)
    print sdbi.get_snp_address("123456_H1234567890-2")
    print sdbi.get_nearest("123456_H1234567890-2")
```
