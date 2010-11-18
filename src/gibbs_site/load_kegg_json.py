#!/usr/bin/env python

import json
import logging

from util import django_utils
from util import kegg

django_utils.SetupDjango()

from gibbs import models


def main():
    parsed_json = json.load(open('kegg_compounds.json'))
    
    for cd in parsed_json:
        try:
            kegg_id = kegg.KeggIdFromInt(cd['CID'])
            c = models.Compound(kegg_id=kegg_id,
                                inchi=cd['InChI'],
                                formula=cd['formula'],
                                mass=cd['mass'])
            # Need to save before setting up many-to-many relationships.
            c.save()

                        
            names = [models.CommonName.GetOrCreate(n)
                     for n in cd['names']]            
            for n in names:
                c.common_names.add(n)
            c.save()
            
        except Exception, msg:
            logging.error(msg)
            continue
        
        
if __name__ == '__main__':
    main()
        
