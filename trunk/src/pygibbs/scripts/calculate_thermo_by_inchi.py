import sys
import csv
from optparse import OptionParser
from pygibbs.groups import GroupContribution, GroupMissingTrainDataError
from toolbox.database import SqliteDatabase
from toolbox.molecule import Molecule
from pygibbs.group_decomposition import GroupDecompositionError
from pygibbs.thermodynamic_constants import default_I, default_pH, default_pMg,\
    default_T

def MakeOpts():
    """Returns an OptionParser object with all the default options."""
    opt_parser = OptionParser()
    opt_parser.add_option("-c", "--csv_input_filename",
                          dest="csv_input_filename",
                          default="../data/thermodynamics/bigg_inchi.csv",
                          help="input CSV file with InChIs")
    opt_parser.add_option("-o", "--csv_output_filename",
                          dest="csv_output_filename",
                          default="../res/bigg_thermo.csv",
                          help="output CSV file with dGs and groups")
    return opt_parser

def CalculateThermo():
    options, _ = MakeOpts().parse_args(sys.argv)

    db = SqliteDatabase('../res/gibbs.sqlite')
    G = GroupContribution(db=db)
    G.init()

    ignore_protonations = False

    print "Reading InChIs from %s, and calculating energies" % options.csv_input_filename
    rowdicts = []
    for row in csv.DictReader(open(options.csv_input_filename, 'r')):
        id = row["ID"]
        if "InChI" in row:
            inchi = row["InChI"]
            diss_table = Molecule._GetDissociationTable(inchi, format='inchi')
        elif "smiles" in row:
            smiles = row["smiles"]
            diss_table = Molecule._GetDissociationTable(smiles, format='smiles')
        else:
            raise Exception("There must be one molecular ID column: InChI or smiles")
        
        try:
            mol = diss_table.GetMostAbundantMol(pH=default_pH, I=default_I, 
                                                pMg=default_pMg, T=default_T)
            
            decomposition = G.Mol2Decomposition(mol, 
                ignore_protonations=ignore_protonations)
            groupvec = decomposition.AsVector()
            dG0 = G.groupvec2val(groupvec)
            nH = decomposition.Hydrogens()
            nMg = decomposition.Magnesiums()
            if nMg > 0:
                raise Exception("Magnesium ions are not allowed here")

            diss_table.SetFormationEnergyByNumHydrogens(dG0, nH, nMg)
            pmap = diss_table.GetPseudoisomerMap()
            
            for p_nH, p_z, p_nMg, p_dG0 in pmap.ToMatrix():
                rowdicts.append({'ID':id, 'nH':p_nH, 'charge':p_z, 'nMg':p_nMg,
                                 'dG0':p_dG0, 'groupvec':str(groupvec),
                                 'error':None})
        except GroupDecompositionError:
            rowdicts.append({'ID':id, 'nH':None, 'charge':None, 'nMg':None,
                             'dG0':None, 'groupvec':None,
                             'error':"cannot decompose"})
        except GroupMissingTrainDataError:
            rowdicts.append({'ID':id, 'nH':None, 'charge':None, 'nMg':None,
                             'dG0':None, 'groupvec':str(groupvec),
                             'error':"missing training data"})
        
    print "writing results to %s" % options.csv_output_filename
    csv_writer = csv.DictWriter(open(options.csv_output_filename, 'w'), 
        ['ID', 'error', 'nH', 'nMg', 'charge', 'dG0', 'groupvec'])
    csv_writer.writeheader()
    csv_writer.writerows(rowdicts)
    
if __name__ == '__main__':
    CalculateThermo()