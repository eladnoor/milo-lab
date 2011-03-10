#!/usr/bin/env python

import json
import logging

from util import django_utils
from util import kegg

django_utils.SetupDjango()

from gibbs import models
        
         
def MirrorCompounds(kegg_json_filename='data/mirror_compounds.json'):
    parsed_json = json.load(open(kegg_json_filename))
    
    for cd in parsed_json:
        to_cid = cd['CID']
        from_cid = cd['mirror from']
        
        to_compound = models.Compound.objects.get(kegg_id=to_cid)
        from_compound = models.Compound.objects.get(kegg_id=from_cid)
        
        # Mirror everything but the KEGG ID
        to_compound.species.clear()
        for specie in from_compound.species.all():
            to_compound.species.add(specie)
            
        to_compound.inchi = from_compound.inchi
        to_compound.formula = from_compound.formula
        to_compound.note = from_compound.note
        to_compound.preferred_name = from_compound.preferred_name
        to_compound.mass = from_compound.mass
        to_compound.num_electrons = from_compound.num_electrons
        to_compound.no_dg_explanation = from_compound.no_dg_explanation
        
        to_compound.save()


if __name__ == '__main__':
    MirrorCompounds()
