DB=$1
echo "Dropping DB: $DB"
dropdb -U ulf -h 158.119.147.72 $DB
echo "Creating DB: $DB"
createdb -U ulf -h 158.119.147.72 $DB
echo "Creating tables in $DB"
psql -U ulf -h 158.119.147.72 $DB < setup_snapper3_db.sql
echo "Finished!"
