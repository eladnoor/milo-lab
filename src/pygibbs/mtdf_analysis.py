#!/usr/bin/python

import logging
import sys

from optparse import OptionParser
from os import path

from metabolic_modelling import bounds
from metabolic_modelling import mtdf_optimizer
from metabolic_modelling import stoich_model
from metabolic_modelling import thermodynamic_data
from pygibbs import kegg
from pygibbs import thermodynamic_estimators
from toolbox import database
from toolbox import util



def MakeOpts():
    """Returns an OptionParser object with all the default options."""
    opt_parser = OptionParser()
    opt_parser.add_option("-k", "--kegg_database_location", 
                          dest="kegg_db_filename",
                          default="../data/public_data.sqlite",
                          help="The KEGG database location")
    opt_parser.add_option("-d", "--database_location", 
                          dest="db_filename",
                          default="../res/gibbs.sqlite",
                          help="The Thermodynamic database location")
    opt_parser.add_option("-s", "--thermodynamics_source",
                          dest="thermodynamics_source",
                          type="choice",
                          choices=thermodynamic_estimators.EstimatorNames(),
                          default="merged",
                          help="The thermodynamic data to use")
    opt_parser.add_option("-i", "--input_filename",
                          dest="input_filename",
                          default="../data/thermodynamics/pathways.txt",
                          help="The file to read for pathways to analyze.")
    opt_parser.add_option("-o", "--output_filename",
                          dest="output_filename",
                          default='../res/thermo_analysis/report.html',
                          help="Where to write output to.")
    return opt_parser

def Main():
    options, _ = MakeOpts().parse_args(sys.argv)
    estimators = thermodynamic_estimators.LoadAllEstimators()
    
    input_filename = path.abspath(options.input_filename)
    if not path.exists(input_filename):
        logging.fatal('Input filename %s doesn\'t exist' % input_filename)
        
    print 'Will read pathway definitions from %s' % input_filename
    
    db_loc = options.db_filename
    print 'Reading from DB %s' % db_loc
    db = database.SqliteDatabase(db_loc)

    thermo = estimators[options.thermodynamics_source]
    print "Using the thermodynamic estimations of: " + thermo.name
    
    # Create a bounds instance
    kegg_instance = kegg.Kegg.getInstance()
    kegg_bounds = kegg_instance.cid2bounds
    my_bounds = bounds.Bounds()
    my_bounds.AddOldStyleBounds(kegg_bounds)

    """
    output_filename = path.abspath(options.output_filename)
    print 'Will write output to %s' % output_filename
    dirname = path.dirname(output_filename)
    if not path.exists(dirname):
        print 'Making output directory %s' % dirname
        util._mkdir(dirname)
    """
    
    print 'Executing MTDF analysis'


if __name__ == "__main__":
    Main()
    
