#!/usr/bin/python

from util import django_utils

django_utils.SetupDjango()

import load_additional_data
import load_formation_energies
import load_kegg_json


def main():
    print 'Loading KEGG data'
    load_kegg_json.LoadAllKeggData()
    
    print 'Loading corrections/additions to KEGG'
    load_additional_data.LoadAdditionalCompoundData()
    
    print 'Loading formation energies'
    load_formation_energies.LoadAllFormationEnergies()
    
    
    
if __name__ == '__main__':
    main()
 