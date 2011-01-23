#!/usr/bin/python

from util import django_utils

django_utils.SetupDjango()

import load_compound_mappings
import load_formation_energies
import load_kegg_json


def main():
    load_kegg_json.LoadKeggJson()
    load_formation_energies.LoadAllFormationEnergies()
    load_compound_mappings.LoadEquivalentCompounds()
 