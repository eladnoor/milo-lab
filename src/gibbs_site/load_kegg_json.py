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
            cid = cd['CID']
            formula = cd['formula']
            mass = float(cd['mass'])
            inchi = cd['InChI']
            num_electrons = cd.get('num_electrons')
            c = models.Compound(kegg_id=cid,
                                formula=formula,
                                inchi=inchi,
                                mass=mass)
            
            if num_electrons != None:
                c.num_electrons = int(num_electrons)
                
            # Need to save before setting up many-to-many relationships.
            c.save()
                        
            names = [models.CommonName.GetOrCreate(n)
                     for n in cd['names']]            
            for n in names:
                c.common_names.add(n)
            c.save()        
        except Exception, e:
            logging.error('Error parsing cid %s', cid)
            logging.debug(e)
            continue
        
if __name__ == '__main__':
    main()
        
