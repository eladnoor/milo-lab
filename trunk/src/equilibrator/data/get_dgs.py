#!/usr/bin/python

import csv
import sys

from util import django_utils
from optparse import OptionParser

django_utils.SetupDjango()

from gibbs import models
from gibbs import constants


def GenFormationEnergyData(pH=constants.DEFAULT_PH,
                          ionic_strength=constants.DEFAULT_IONIC_STRENGTH):
    """Returns a list of dictionaries of compound formation energies.
    
    TODO(flamholz): return data from multiple sources per compound when possible.
    
    Args:
        pH: the pH.
        ionic_strength: the ionic strength.
    """
    dicts = {}
    for compound in models.Compound.objects.all():
        dG = compound.DeltaG(pH=pH, ionic_strength=ionic_strength)
        d = {'KEGG ID': compound.kegg_id,
             'Name': compound.FirstName(),
             'dG (kJ/mol)': dG,
             'InChI': compound.inchi,
             'Source': compound.dg_source}
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
    opt_parser.add_option("-o", "--output_filename", dest="output_filename",
                          type="string", default="formation_energies.csv",
                          help="The name of the file to write output to.")
    return opt_parser


def main():
    options, _ = MakeOpts().parse_args(sys.argv)
    print 'Using pH = %.2f and ionic strength = %.3f' % (options.pH,
                                                         options.ionic_strength)
    print 'Will write output to %s' % options.output_filename

    dicts = GenFormationEnergyData(pH=options.pH,
                                   ionic_strength=options.ionic_strength)
    csv_data = sorted(dicts, key=lambda x: (x['KEGG ID'], x['Source']))
    csv_file = open(options.output_filename, 'w')
    writer = csv.DictWriter(csv_file, ['KEGG ID', 'Name', 'dG (kJ/mol)', 'InChI', 'Source'])
    writer.writerows(csv_data)
        
    csv_file.close()
    
    print 'Done.'
            

if __name__ == '__main__':
    main()
                
                
