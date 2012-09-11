#!/usr/bin/python

import logging
import numpy as np
import sys

from os import path

from pygibbs.metabolic_modelling import kinetic_data
from pygibbs.metabolic_modelling import feasible_concentrations_iterator
from pygibbs.metabolic_modelling import protein_optimizer
from pygibbs.metabolic_modelling import thermodynamic_data
from pygibbs import kegg
from pygibbs import thermodynamic_estimators
from pygibbs import pathway
from pygibbs import templates
from toolbox import util
from toolbox import stats
from argparse import ArgumentParser


def MakeOpts():
    """Returns an OptionParser object with all the default args."""
    parser = ArgumentParser()
    parser.add_argument("-k", "--kegg_database_location", 
                          dest="kegg_db_filename",
                          default="../data/public_data.sqlite",
                          help="The KEGG database location")
    parser.add_argument("-d", "--database_location", 
                          dest="db_filename",
                          default="../res/gibbs.sqlite",
                          help="The Thermodynamic database location")
    parser.add_argument("-s", "--thermodynamics_source",
                          dest="thermodynamics_source",
                          type="choice",
                          choices=thermodynamic_estimators.EstimatorNames(),
                          default="merged",
                          help="The thermodynamic data to use")
    parser.add_argument("-n", "--kinetics_filename",
                          dest="kinetics_filename",
                          default=None,
                          help="The kinetics data file to use or None.")
    parser.add_argument("-i", "--input_filename",
                          dest="input_filename",
                          default="../data/thermodynamics/pathways.txt",
                          help="The file to read for pathways to analyze.")
    parser.add_argument("-o", "--output_dir",
                          dest="output_dir",
                          default='../res/protein_analysis/',
                          help="Where to write output to.")
    return parser


def Main():
    np.seterr('raise')
    parser = MakeOpts()
    args = parser.parse_args()
    estimators = thermodynamic_estimators.LoadAllEstimators()
    
    input_filename = path.abspath(args.input_filename)
    if not path.exists(input_filename):
        logging.fatal('Input filename %s doesn\'t exist' % input_filename)
        
    print 'Will read pathway definitions from %s' % input_filename

    # Make thermodynamic and kinetic data containers
    thermo = estimators[args.thermodynamics_source]
    print "Using the thermodynamic estimations of: " + thermo.name
    thermo_data = thermodynamic_data.WrapperThermoData(thermo)
    
    # Fetch kinetic data.
    kin_data = kinetic_data.UniformKineticData(kcat=200, km=2e-4, mass=40)
    if args.kinetics_filename is not None:
        print 'Parsing kinetic data from', args.kinetics_filename
        kin_data = kinetic_data.KineticDataWithDefault.FromArrenFile(
            args.kinetics_filename)
        
    """
    kin_data = kinetic_data.KineticDataWithDefault.FromFiles(
        '../data/enzymatics/glycolytic_pathway_enzymes_kcat.csv',
        '../data/enzymatics/glycolytic_pathway_enzymes_km.csv')
    kin_data.SetDefaultKcat(100)
    kin_data.SetDefaultKM(1e-4)
    kin_data.SetDefaultMass(35)
    """
    
    # Create a kegg instance
    kegg_instance = kegg.Kegg.getInstance()

    # Create output directories
    out_dir = args.output_dir
    if not path.exists(out_dir):
        util._mkdir(out_dir)
    pathgraph_dir = path.join(out_dir, 'pathway_graphs/')
    util._mkdir(pathgraph_dir)
    
    print 'Executing Protein Cost analysis'
    pathway_iterator = pathway.KeggPathwayIterator.FromFilename(input_filename)
    results = []
    for pathway_data in pathway_iterator:
        if pathway_data.skip:
            print 'Skipping pathway', pathway_data.name
            continue
        
        print 'Analyzing pathway', pathway_data.name
                
        model = pathway_data.GetStoichiometricModel(kegg_instance)
        model_bounds = pathway_data.GetBounds()

        opt = protein_optimizer.ProteinOptimizer(model, thermo_data, kin_data)
        it = feasible_concentrations_iterator.FeasibleConcentrationsIterator(
            model, thermo_data, model_bounds)
        
        # Now solve with the default initial conditions.
        success = None
        result = None
        optima = []
        for i, x0 in enumerate(it):
            result = opt.FindOptimum(model_bounds, initial_concentrations=x0)
            status = result.status
            print '\t%s optimization %d' % (pathway_data.name, i)
            if status.failure:          
                print '\tFailed to optimize', pathway_data.name
                print '\t%s' % status
            elif status.infeasible:      
                print '\t', pathway_data.name, 'is infeasible!'
                print '\t%s' % status
            else:
                print '\t*Optimization successful'
                optima.append(result.opt_val)
                if not success:
                    success = result
                elif result.opt_val < success.opt_val:
                    success = result
        
        mean, error = None, None
        if optima:
            try:
                mean, error = stats.MeanWithConfidenceInterval(optima)
            except Exception, e:
                mean, error = None, None
                print optima
        result_dict = {'result': None,
                       'num_optima': len(optima),
                       'mean_opt': mean,
                       'error': error}
        
        if success is not None:
            success.WriteAllGraphs(pathgraph_dir)
            result_dict['result'] = success
        
            cost = success.opt_val
            if cost is not None:
                print '\t*Protein Cost for', pathway_data.name, '= %.2g' % cost
            if optima:
                print 'Found', len(optima), 'near-optima for', pathway_data.name 
                optima = np.array(optima)
                mean_opt = np.mean(optima)
                mean_diff = np.mean(np.abs(optima - mean_opt))
                print 'Mean optimum', mean_opt
                print 'Mean diff from mean', mean_diff
                print 'Percent diff %s%%' % (100*mean_diff / mean_opt)
                print 'StdDev opt', np.std(optima)
        else:
            # Use default conditions to show the failure
            res = opt.FindOptimum(model_bounds)
            result_dict['result'] = res            
        results.append(result_dict)
    
    output_filename = path.join(out_dir, 'results.html')
    print 'Writing output to', output_filename
    template_data = {'analysis_type': 'Protein Cost',
                     'kinetic_data': kin_data,
                     'results':results}
    templates.render_to_file('protein_optimization_results.html',
                             template_data,
                             output_filename)
    

if __name__ == "__main__":
    Main()
    
