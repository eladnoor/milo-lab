#/bin/bash

echo "You should run this script only after making changes to data/kegg/kegg_additions.csv."
echo "This script will run all necessary procedures to update the local eQuilibrator database."
echo "WARNING: It must not be interrupted!!!"
echo -n "Are you sure you wish to proceed (Y/N) ?"
read text
if [ $text == 'Y' ]; then
    python pygibbs/kegg.py
    python pygibbs/dissociation_constants.py -k
    python pygibbs/groups.py -e
    python pygibbs/scripts/export_kegg_data.py
    gunzip -c ../res/kegg_compounds.json.gz > equilibrator/data/kegg_compounds.json
    cd equilibrator
    python load_database.py
    cd ..
else
    echo "Okay, see you next time"
fi
