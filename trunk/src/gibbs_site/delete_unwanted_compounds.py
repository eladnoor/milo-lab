#!/usr/bin/env python

import json
import logging

from util import django_utils

django_utils.SetupDjango()

from gibbs import models


def DeleteUnwantedCompounds(json_filename='data/compounds_to_delete.json'):
    parsed_json = json.load(open(json_filename))

    for cid in parsed_json:
        try:
            compound = models.Compound.objects.get(kegg_id=cid)
            compound.delete()
        except Exception, e:
            logging.error('Error deleting cid %s', cid)
            continue

    
if __name__ == '__main__':
    DeleteUnwantedCompounds()