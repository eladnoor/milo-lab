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
from pygibbs import pathway
from pygibbs import templates
from toolbox import database

import numpy as np


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
    thermo_data = thermodynamic_data.WrapperThermoData(thermo)
    
    # Create a bounds instance
    kegg_instance = kegg.Kegg.getInstance()

    
    print 'Executing MTDF analysis'
    pathway_iterator = pathway.KeggPathwayIterator.FromFilename(input_filename)
    mtdf_results = []
    for pathway_data in pathway_iterator:
        if pathway_data.skip:
            print 'Skipping pathway', pathway_data.name
            continue
        
        
        print 'Analyzing pathway', pathway_data.name
        
        model = pathway_data.GetStoichiometricModel(kegg_instance)
        model_bounds = pathway_data.GetBounds()
        
        mtdf_opt = mtdf_optimizer.MTDFOptimizer(model, thermo_data)
        mtdf_result = mtdf_opt.FindMTDF(model_bounds)
        mtdf_results.append(mtdf_result)
        
        mtdf = mtdf_result.mtdf
        print '\tMTDF for', pathway_data.name, '= %.2g' % mtdf
    
    
    output_filename = path.abspath(options.output_filename)
    template_data = {'mtdf_results':mtdf_results}
    templates.render_to_file('mtdf_results.html',
                             template_data,
                             output_filename)
    

if __name__ == "__main__":
    Main()
    
