#!/usr/bin/python

import csv
import logging

from util import django_utils
from util import kegg

django_utils.SetupDjango()

from gibbs import models


def main():    
    reader = csv.DictReader(open('formation.csv'))
    
    species = []
    values = []
    for i, line in enumerate(reader):
        kegg_id = kegg.KeggIdFromInt(int(line['cid']))

        estimated = line['estimated'] == '1'
        source = models.ValueSource.Alberty()
        if estimated:
            source = models.ValueSource.GroupContribution()

        specie = models.Specie(kegg_id=kegg_id,
                               number_of_hydrogens=line['nH'],
                               net_charge=line['z'],
                               formation_energy=line['dG0'],
                               formation_energy_source=source)
        specie.save()
        
        c = models.Compound.objects.get(kegg_id=kegg_id)
        c.species.add(specie)
        c.save()            


if __name__ == '__main__':
    main()