#!/usr/bin/env python

import json
import logging

from util import django_utils

django_utils.SetupDjango()

from gibbs import models


def LoadEquivalentCompounds(
    replacements_json_filename='data/compound_replacement_mapping.json'):
    """Loads mappings between compounds."""
    parsed_json = json.load(open(replacements_json_filename))
    
    for mapping in parsed_json:
        from_kegg_id = mapping['FROM ID']
        to_kegg_id = mapping['TO ID']
        
        try:
            from_compound = models.Compound.objects.get(kegg_id=from_kegg_id)
            to_compound = models.Compound.objects.get(kegg_id=to_kegg_id)
            from_compound.replace_with = to_compound
            from_compound.save()
        except Exception, e:
            logging.error(e)
            continue


if __name__ == '__main__':
    LoadEquivalentCompounds()