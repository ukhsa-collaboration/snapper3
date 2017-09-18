# snapper3

A beta version of the software is provided in **/phengs/hpc_software/phe/snapper3/3-0**. To use it:

```bash
module load phe/snapper3/3-0
```

Prerequisites/dependencies:
- Python >= 2.7.6
- psycopg2 >= 2.5.2

The postgres server must be running postgres >=9.6 and the contrib package must be installed so that the intarray extension can be created for each database.

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

To create a new database, use the reset.sh script and the use the add_reference subcommand of snapper3.py

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





To add a new sample to the database
