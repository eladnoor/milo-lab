#!/usr/bin/env python

import json
import logging

from util import django_utils

django_utils.SetupDjango()

from gibbs import models

DEFAULT_CITATION_DATA_FILENAME = 'data/citation_data.json'


def CheckData(filenames=(DEFAULT_CITATION_DATA_FILENAME,)):
    for json_fname in filenames:
        json.load(open(json_fname))


def LoadCitationData(json_filename=DEFAULT_CITATION_DATA_FILENAME):
    parsed_json = json.load(open(json_filename))

    for cd in parsed_json:
        try:
            name = cd['name']
            ref = cd['ref']
            
            source = models.ValueSource(name=name, citation=ref)
            source.save()
        except Exception, e:
            logging.error('Error parsing reference %s', cd)
            logging.error(e)
            continue


if __name__ == '__main__':
    LoadCitationData()