#!/usr/bin/python

import csv
import sys
import json

from util import django_utils
from optparse import OptionParser

django_utils.SetupDjango()

from gibbs import models
from gibbs import constants

# Column names
KEGG_ID = '!MiriamID::urn:miriam:kegg.compound'
NAME = '!Name'
INCHI = '!InChI'
SOURCE = '!Source'
FORMATION = '!FormationEnergy'
ROW_ORDER = [NAME, KEGG_ID, INCHI, FORMATION, SOURCE]


def GenFormationEnergyData(pH=constants.DEFAULT_PH,
                          ionic_strength=constants.DEFAULT_IONIC_STRENGTH):
    """Returns a list of dictionaries of compound formation energies.
    
    TODO(flamholz): return data from multiple sources per compound when possible.
    
    Args:
        pH: the pH.
        ionic_strength: the ionic strength.
    """
    dicts = []
    for compound in models.Compound.objects.all():
        dG = compound.DeltaG(pH=pH, ionic_strength=ionic_strength)
        if dG:
            dG = round(dG, 3)
        d = {KEGG_ID: compound.kegg_id,
             NAME: compound.FirstName(),
             FORMATION: dG,
             INCHI: compound.inchi,
             SOURCE: compound.dg_source.name}
        dicts.append(d)
    return dicts


def MakeOpts():
    """Returns an OptionParser object with all the default options."""
    opt_parser = OptionParser()
    opt_parser.add_option("-p", "--ph", dest="pH", type="float",
                          default=constants.DEFAULT_PH, help="The pH")
    opt_parser.add_option("-i", "--ionic_strength", dest="ionic_strength",
                          type="float", default=constants.DEFAULT_IONIC_STRENGTH,
                          help="The ionic strength, M")
    opt_parser.add_option("-c", "--output_csv", dest="output_csv",
                          type="string", default="formation_energies.csv",
                          help="The name of the file to write csv output to.")
    opt_parser.add_option("-j", "--output_json", dest="output_json",
                          type="string", default="formation_energies.json",
                          help="The name of the file to write json output to.")
    return opt_parser


def main():
    options, _ = MakeOpts().parse_args(sys.argv)
    print 'Using pH = %.2f and ionic strength = %.3f' % (options.pH,
                                                         options.ionic_strength)
    print 'Will write csv output to %s' % options.output_csv
    print 'Will write json output to %s' % options.output_json

    dicts = GenFormationEnergyData(pH=options.pH,
                                   ionic_strength=options.ionic_strength)
    sorted_data = sorted(dicts, key=lambda x: (x[KEGG_ID], x[SOURCE]))
    csv_file = open(options.output_csv, 'w')
    writer = csv.DictWriter(csv_file, ROW_ORDER)
    writer.writeheader()
    writer.writerows(sorted_data)
    csv_file.close()
    
    json_file = open(options.output_json, 'w')
    json.dump(sorted_data, json_file, sort_keys=True, indent=4)
    json_file.close()
    
    print 'Done.'
            

if __name__ == '__main__':
    main()
                
                
