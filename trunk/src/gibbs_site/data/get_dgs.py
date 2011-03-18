#!/usr/bin/python

from util import django_utils
import csv
import logging

# NOTE(flamholz): This is crappy. We're using the real database for
# a unit test. I wish I knew of a better way.
django_utils.SetupDjango()

from gibbs import models

def main():
    compounds = models.Compound.objects.all()
    
    
    csv_file = open('all_dgs.csv', 'w')
    writer = csv.DictWriter(csv_file, ['KEGG ID', 'Name', 'dG (kJ/mol)', 'source'])
    for compound in compounds:
        dG = compound.DeltaG()
        if dG is None:
            continue
        
        d = {'KEGG ID': compound.kegg_id,
             'Name': compound.FirstName(),
             'dG (kJ/mol)': dG,
             'source': compound.dg_source}
        writer.writerow(d)
        
    csv_file.close()
            

if __name__ == '__main__':
    main()
                
                
