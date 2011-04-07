#!/usr/bin/python

import json
import logging

from util import django_utils
from util import inchi
from util import kegg

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


def GetSource(source_string):
    if not source_string:
        return None
    
    lsource = source_string.lower()
    if lsource.startswith('alberty'):
        return models.ValueSource.Alberty()
    elif lsource.startswith('thauer'):
        return models.ValueSource.Thauer()
    elif lsource.startswith('group'):
        return models.ValueSource.GroupContribution()
    return None


def LoadFormationEnergies(json):
    for cdict in json:
        cid = cdict.get('cid')
        if not cid:
            logging.error('No kegg ID!')
            logging.error(cdict)
            continue
        kegg_id = kegg.KeggIdFromInt(cid)

        try:
            compound = models.Compound.objects.get(kegg_id=kegg_id)
        except Exception, e:
            logging.error(e)
            logging.error(cdict)
            continue
            
        species = cdict.get('species')
        if not species:
            error = cdict.get('error')
            if error:
                compound.no_dg_explanation = error
                compound.save()
                continue
            
        # Override the passed-in source if one is specified in JSON.
        source_string = cdict.get('source')
        my_source = GetSource(source_string)
        if not my_source:
            logging.error('Couldn\'t get source for %s' % my_source)
            logging.error(cdict)
            continue
        
        AddAllSpeciesToCompound(compound, species, my_source)
        compound.save()
            

DEFAULT_DATA_FILE = 'data/pseudoisomers.json'

def CheckData(json_filename=DEFAULT_DATA_FILE):
    json.load(open(json_filename))
    

def LoadAllFormationEnergies(json_filename=DEFAULT_DATA_FILE):
    json_data = json.load(open(json_filename))
    
    print 'Writing Pseudoisomers data'
    LoadFormationEnergies(json_data)


if __name__ == '__main__':
    LoadAllFormationEnergies()