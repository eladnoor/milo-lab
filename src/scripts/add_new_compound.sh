#/bin/bash

echo "You should run this script only after making changes to data/kegg/kegg_additions.csv."
echo "This script will run all necessary procedures to update the local eQuilibrator database."
echo "WARNING: It must not be interrupted!!!"
echo -n "Are you sure you wish to proceed (Y/N)? "
read text
if [ $text == 'Y' ]; then
    python pygibbs/kegg.py
    python pygibbs/dissociation_constants.py -k
    python pygibbs/unified_group_contribution.py -g -o -v -m -t
    python pygibbs/scripts/export_kegg_data.py -s UGC
    pushd equilibrator
    cp -f ../../res/kegg_compounds.json.gz data/kegg_compounds.json.gz
    python manage.py flush --noinput
    python load_database.py
    python manage.py runserver
    popd
else
    echo "Okay, see you next time"
fi
