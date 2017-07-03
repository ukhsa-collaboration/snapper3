HOST=$1
DB=$2
echo "Dropping DB: $DB"
dropdb -U ulf -h $HOST $DB
echo "Creating DB: $DB"
createdb -U ulf -h $HOST $DB
echo "Creating tables in $DB"
psql -U ulf -h $HOST $DB < setup_snapper3_db.sql
echo "Finished!"
