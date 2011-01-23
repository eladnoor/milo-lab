#!/usr/bin/python

import json
import logging

from util import django_utils
from util import inchi

django_utils.SetupDjango()

from gibbs import models


def AddAllSpeciesToCompound(compound, species_dicts, source):
    print 'Writing data from source %s for compound %s' % (source.name,
                                                           compound.kegg_id)
    
    compound.species.clear()
    for sdict in species_dicts:
        specie = models.Specie(kegg_id=compound.kegg_id,
                               number_of_hydrogens=sdict['nH'],
                               number_of_mgs=sdict['nMg'],
                               net_charge=sdict['z'],
                               formation_energy=sdict['dG0_f'],
                               formation_energy_source=source)
        specie.save()
        compound.species.add(specie)
    compound.save()


def LoadFormationEnergies(json, source):
    for cdict in json:
        inchi_str = cdict['inchi']
        if not inchi_str:
            logging.error('No inchi!')
            logging.error(cdict)
            continue
        
        compounds = models.Compound.objects.filter(inchi__exact=inchi_str)
        for c in compounds:
            if source == models.ValueSource.Alberty() or not len(c.species.all()):
                AddAllSpeciesToCompound(c, cdict['species'], source)
            

def LoadAllFormationEnergies(alberty_json_filename='data/alberty.json',
                             gc_json_filename='data/group_contribution.json'):
    alberty_json = json.load(open(alberty_json_filename))
    gc_json = json.load(open(gc_json_filename))
    
    print 'Writing Alberty data'
    LoadFormationEnergies(alberty_json, models.ValueSource.Alberty())
    
    print 'Writing Group Contribution data'
    LoadFormationEnergies(gc_json, models.ValueSource.GroupContribution()) 


if __name__ == '__main__':
    main()