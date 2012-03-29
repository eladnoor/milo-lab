import sys
import csv
from optparse import OptionParser
from toolbox.database import SqliteDatabase
from toolbox.molecule import Molecule, OpenBabelError
from pygibbs.group_decomposition import GroupDecompositionError
import pybel
from pygibbs.unified_group_contribution import UnifiedGroupContribution,\
    UnknownReactionEnergyError
import numpy as np
from pygibbs.dissociation_constants import MissingDissociationConstantError

def MakeOpts():
    """Returns an OptionParser object with all the default options."""
    opt_parser = OptionParser()
    opt_parser.add_option("-i", "--sdf_input_filename",
                          dest="sdf_input_filename",
                          default="../data/metabolic_models/iAF1260.sdf",
                          help="input SDF file with MOL descriptions of multiple compounds")
    opt_parser.add_option("-o", "--csv_output_filename",
                          dest="csv_output_filename",
                          default="../res/iAF1260.csv",
                          help="output CSV file with chemical dG0")
    return opt_parser

def CalculateThermo():
    options, _ = MakeOpts().parse_args(sys.argv)

    if options.csv_output_filename is not None:
        out_fp = open(options.csv_output_filename, 'w')
        print "writing results to %s ... " % options.csv_output_filename
    else:
        out_fp = sys.stdout
    csv_writer = csv.writer(out_fp)
    csv_writer.writerow(['ID', 'error', 'nH', 'nMg', 'charge', 'dG0', 'kernel']) 

    db = SqliteDatabase('../res/gibbs.sqlite', 'w')
    ugc = UnifiedGroupContribution(db)
    ugc.LoadGroups(True)
    ugc.LoadObservations(True)
    ugc.LoadGroupVectors(True)
    ugc.LoadData(True)
    
    result_dict = ugc._GetContributionData(ugc.S.copy(), ugc.cids,
                                           ugc.b.copy(), ugc.anchored)
    
    g_pgc = result_dict['group_contributions']
    P_L_pgc = result_dict['pgc_conservations']

    sdfile = pybel.readfile("sdf", options.sdf_input_filename)
    for m in sdfile:
        try:
            try:
                mol = Molecule.FromOBMol(m.OBMol)
            except OpenBabelError:
                raise UnknownReactionEnergyError("Cannot convert to OBMol object")
            
            mol.title = m.title
            mol.RemoveHydrogens()
            if mol.GetNumAtoms() > 200:
                raise UnknownReactionEnergyError("Compound contains more than 200 atoms (n = %d)" % mol.GetNumAtoms())
            
            try:
                decomposition = ugc.group_decomposer.Decompose(mol, 
                                        ignore_protonations=False, strict=True)
            except GroupDecompositionError:
                raise UnknownReactionEnergyError("cannot decompose")
            
            groupvec = decomposition.AsVector()
            gv = np.matrix(groupvec.Flatten())
            dG0 = float(g_pgc * gv.T)
            nH = decomposition.Hydrogens()
            nMg = decomposition.Magnesiums()
            ker = list((P_L_pgc * gv.T).round(10).flat)
            try:
                diss_table = mol.GetDissociationTable()
                diss_table.SetFormationEnergyByNumHydrogens(
                        dG0=dG0, nH=nH, nMg=nMg)
            except MissingDissociationConstantError:
                raise UnknownReactionEnergyError("missing pKa data")
            pmap = diss_table.GetPseudoisomerMap()
            for p_nH, p_z, p_nMg, p_dG0 in pmap.ToMatrix():
                csv_writer.writerow([m.title, None, p_nH, p_z, p_nMg, round(p_dG0, 1), str(ker)])

        except UnknownReactionEnergyError as e:
            csv_writer.writerow([m.title, str(e), None, None, None, None, None])
        
        out_fp.flush()
if __name__ == '__main__':
    CalculateThermo()