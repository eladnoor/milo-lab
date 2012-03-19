#!/usr/bin/python

import logging
import sys
import numpy as np
import pylab

from optparse import OptionParser
from os import path
from matplotlib.font_manager import FontProperties

from pygibbs.metabolic_modelling import kinetic_data
from pygibbs.metabolic_modelling import protein_optimizer
from pygibbs.metabolic_modelling import mtdf_optimizer
from pygibbs.metabolic_modelling import thermodynamic_data
from pygibbs.thermodynamic_constants import default_RT as RT
from pygibbs import kegg
from pygibbs import thermodynamic_estimators
from pygibbs import pathway
from toolbox import util

LEGEND_FONT = FontProperties(size=8)

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
    opt_parser.add_option("-o", "--output_dir",
                          dest="output_dir",
                          default='../res/protein_analysis/',
                          help="Where to write output to.")
    return opt_parser


def Main():
    options, _ = MakeOpts().parse_args(sys.argv)
    estimators = thermodynamic_estimators.LoadAllEstimators()
    
    input_filename = path.abspath(options.input_filename)
    if not path.exists(input_filename):
        logging.fatal('Input filename %s doesn\'t exist' % input_filename)
        
    print 'Will read pathway definitions from %s' % input_filename

    # Make thermodynamic and kinetic data containers
    thermo = estimators[options.thermodynamics_source]
    print "Using the thermodynamic estimations of: " + thermo.name
    thermo_data = thermodynamic_data.WrapperThermoData(thermo)
    
    # Uniform kinetic data
    kin_data = kinetic_data.UniformKineticData(kcat=100, km=1e-4)
    
    # Create a kegg instance
    kegg_instance = kegg.Kegg.getInstance()

    # Create output directories
    out_dir = options.output_dir
    if not path.exists(out_dir):
        util._mkdir(out_dir)
    pathgraph_dir = path.join(out_dir, 'pathway_graphs/')
    util._mkdir(pathgraph_dir)
    
    print 'Executing Protein Cost analysis'
    pathway_iterator = pathway.KeggPathwayIterator.FromFilename(input_filename)
    mtdfs = []
    protein_scores = []
    names = []
    num_atp = []
    path_lengths = []
    for pathway_data in pathway_iterator:
        if pathway_data.skip:
            print 'Skipping pathway', pathway_data.name
            continue
        
        print 'Analyzing pathway', pathway_data.name
                
        model = pathway_data.GetStoichiometricModel(kegg_instance)
        model_bounds = pathway_data.GetBounds()
        
        protein_opt = protein_optimizer.ProteinOptimizer(model, thermo_data, kin_data)
        mtdf_opt = mtdf_optimizer.MTDFOptimizer(model, thermo_data)
        
        # Solve MTDF.
        mtdf_res = mtdf_opt.FindMTDF(model_bounds)
        mtdf_status = mtdf_res.status
        if mtdf_status.IsFailure() or mtdf_status.IsInfeasible():
            print '\tFailed to optimize', pathway_data.name
            continue

        # Solve protein.
        protein_res = protein_opt.FindOptimum(model_bounds)
        protein_status = protein_res.status
        if protein_status.IsFailure() or protein_status.IsInfeasible():          
            print '\tFailed to optimize', pathway_data.name
            continue

        mtdfs.append(mtdf_res.opt_val)
        protein_scores.append(protein_res.opt_val)
        names.append(model.name)
        
        net_reaction = mtdf_res.net_reaction.sparse
        atp_produced = net_reaction.get(2, 0)
        num_atp.append(atp_produced)
        path_lengths.append(len(mtdf_res.reaction_ids))
        
        pylab.figure()
        pylab.title(model.name)
        dGr0_tag = mtdf_res.dGr0_tag.flatten().tolist()
        dgmtdf = mtdf_res.dGr_tag.flatten().tolist()
        dgprotein = protein_res.dGr_tag.flatten().tolist()
        dgbio = mtdf_res.dGr_bio.flatten().tolist()
        dg0_profile = np.cumsum([0] + dGr0_tag)
        dgmtdf_profile = np.cumsum([0] + dgmtdf)
        dgprotein_profile = np.cumsum([0] + dgprotein)
        dgbio_profile = np.cumsum([0] + dgbio)
        
        rxn_range = pylab.arange(len(mtdf_res.reaction_ids) + 1)
        pylab.plot(rxn_range, dg0_profile, 'b--',
                   linewidth=2, label='Standard Conditions')
        pylab.plot(rxn_range, dgbio_profile, 'c--',
                   linewidth=2, label='Biological Conditions')
        mtdf_label = 'MTDF Optimized (MTDF = %.2g kJ/mol)' % mtdf_res.opt_val
        pylab.plot(rxn_range, dgmtdf_profile, 'r-',
                   linewidth=2, label=mtdf_label)
        pc_label = 'Protein Optimized (Cost = %.2g)' % protein_res.opt_val
        pylab.plot(rxn_range, dgprotein_profile, 'g-',
                   linewidth=2, label=pc_label)
        pylab.xticks(rxn_range[:-1] + 0.5, mtdf_res.reaction_ids)
        pylab.xlabel('Reaction step')
        pylab.ylabel('Cumulative dG (kJ/mol)')
        pylab.legend(loc='upper right', prop=LEGEND_FONT)
    
    pylab.figure()
    pylab.plot(num_atp, protein_scores, 'b.')
    #pylab.xlabel('MTDF (kJ/mol)')
    pylab.xlabel('Net ATP Production')
    pylab.ylabel('Protein Cost')
    for x,y,s in zip(num_atp, protein_scores, names):
        pylab.text(x, y, s, fontsize=10)
    
    max_protein = np.max(protein_scores)
    pylab.plot([0,0], [0,max_protein], 'r--', label='0 ATP Produced')
    pylab.plot([1,1], [0,max_protein], 'g--', label='1 ATP Produced')
    pylab.plot([2,2], [0,max_protein], 'b--', label='2 ATP Produced')
    
    #pylab.yscale('log')
    pylab.xticks([])
    pylab.xlim((-1, 3))
    pylab.legend()
    
    odbs = np.tanh(np.array(mtdfs) / (2*RT))
    
    pylab.figure()
    pylab.plot(protein_scores, odbs, 'b.')
    pylab.xlabel('Protein Cost')
    pylab.ylabel('ODB (unitless)')
    
    #for x,y,s in zip(protein_scores, length_scaled_cost, names):
    #    pylab.text(x, y, s, fontsize=10)  
    pylab.show()
    

    

if __name__ == "__main__":
    Main()
    
