#!/usr/bin/env python

import json
import logging

from util import django_utils

django_utils.SetupDjango()

from gibbs import models


def GetOrCreateNames(names_list):
    """Find all the names in the database.
    
    Create them if they are not present.
    """
    return [models.CommonName.GetOrCreate(n)
            for n in names_list]


def LoadAdditionalCompoundData(json_filename='data/additional_compound_data.json'):
    parsed_json = json.load(open(json_filename))

    for cd in parsed_json:
        try:
            cid = cd['CID']
            compound = models.Compound.objects.get(kegg_id=cid)
            
            note = cd.get('note')
            preferred_name = cd.get('preferred name')
            if note:
                compound.note = note
            if preferred_name:
                compound.preferred_name = preferred_name
            
            names = cd.get('names')
            if names:
                for n in GetOrCreateNames(names):
                    compound.common_names.add(n)
                    
            compound.save()
        except Exception, e:
            logging.error('Error parsing cid %s', cid)
            logging.error(e)
            continue


def LoadAllAdditionalData():
    LoadAdditionalCompoundData()

    
if __name__ == '__main__':
    LoadAdditionalCompoundData()