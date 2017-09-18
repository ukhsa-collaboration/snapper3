USER=$1
HOST=$2
DB=$3
echo "Dropping DB: $DB"
dropdb -U $USER -h $HOST $DB
echo "Creating DB: $DB"
createdb -U $USER -h $HOST $DB
echo "Creating tables in $DB"
psql -U $USER -h $HOST $DB < setup_snapper3_db.sql
echo "Creating functions in $DB"
psql -U $USER -h $HOST $DB < add_psql_functions.sql
echo "Finished!"
