#!/usr/bin/python

from util import django_utils

django_utils.SetupDjango()

import load_additional_data
import load_citation_data
import load_kegg_json


def main():
    load_citation_data.CheckData()
    load_kegg_json.CheckData()
    load_additional_data.CheckData()

    print 'Loading citation data'
    load_citation_data.LoadCitationData()
        
    print 'Loading KEGG data'
    load_kegg_json.LoadAllKeggData()
    
    print 'Loading corrections/additions to KEGG'
    load_additional_data.LoadAdditionalCompoundData()    
    
    
if __name__ == '__main__':
    main()
 