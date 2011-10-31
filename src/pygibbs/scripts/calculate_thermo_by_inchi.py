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
    opt_parser.add_option("-s", "--smiles",
                          dest="smiles",
                          default=None,
                          help="a single molecule input in SMILES format")
    opt_parser.add_option("-c", "--csv_input_filename",
                          dest="csv_input_filename",
                          default=None,
                          help="input CSV file with InChIs")
    opt_parser.add_option("-o", "--csv_output_filename",
                          dest="csv_output_filename",
                          default=None,
                          help="output CSV file with dGs and groups")
    opt_parser.add_option("-b", "--biochemical",
                          dest="biochemical",
                          default=False,
                          action="store_true",
                          help="calculate the biochemical (transformed) energy at standard conditions")
    opt_parser.add_option("-p", "--ph", action="store", type="float",
                          dest="pH", default=default_pH,
                          help="The width of the range of pH values to plot")
    opt_parser.add_option("-i", "--ionic_strength", action="store", type="float",
                          dest="I", default=default_I,
                          help="The Ionic strength of the solution (in M)")
    opt_parser.add_option("-m", "--pmg", action="store", type="float",
                          dest="pMg", default=default_pMg,
                          help="The pMg of the solution")
    opt_parser.add_option("-t", "--temperature", action="store", type="float",
                          dest="T", default=default_T,
                          help="The T of the solution (in K)")
    return opt_parser

def CalculateThermo():
    parser = MakeOpts()
    options, _ = parser.parse_args(sys.argv)
    pH, I, pMg, T = options.pH, options.I, options.pMg, options.T

    db = SqliteDatabase('../res/gibbs.sqlite')
    G = GroupContribution(db=db)
    G.init()
    ignore_protonations = False

    list_of_mols = []
    if options.smiles:
        list_of_mols.append({'id':options.smiles, 'mol':options.smiles,
            'format':'smiles'})
    elif options.csv_input_filename:
        for row in csv.DictReader(open(options.csv_input_filename, 'r')):
            if "InChI" in row:
                list_of_mols.append({'id':row["ID"], 'mol':row["InChI"],
                                     'format':'inchi'})
            elif "smiles" in row:
                list_of_mols.append({'id':row["ID"], 'mol':row["smiles"],
                                     'format':'smiles'})
            else:
                raise Exception("There must be one molecular ID column: InChI or smiles")
    else:
        parser.error("must use either -s or -c option")
    
    if options.biochemical:
        print ("Calculating biochemical formation energies for %s compounds" 
               " at pH = %.1f, I = %.2f, pMg = %.1f, T = %.2f" %  
               (len(list_of_mols), pH, I, pMg, T))
    else:
        print ("Calculating chemical formation energies for %s compounds" % 
               len(list_of_mols))
    
    rowdicts = []
    for mol_dict in list_of_mols:
        mol_id = mol_dict['id']
        diss_table = Molecule._GetDissociationTable(mol_dict['mol'],
                                                    format=mol_dict['format'])
        try:
            mol = diss_table.GetMostAbundantMol(pH, I, pMg, T)
            decomposition = G.Mol2Decomposition(mol, 
                ignore_protonations=ignore_protonations)
            groupvec = decomposition.AsVector()
            dG0 = G.groupvec2val(groupvec)
            nH = decomposition.Hydrogens()
            nMg = decomposition.Magnesiums()
            diss_table.SetFormationEnergyByNumHydrogens(dG0, nH, nMg)
            pmap = diss_table.GetPseudoisomerMap()
            
            if options.biochemical:
                dG0_prime = pmap.Transform(pH, pMg, I, T)
                rowdicts.append({'ID':mol_id, 'pH':pH, 'I':I, 'pMg':pMg,
                                 'dG0\'':"%.1f" % dG0_prime, 'groupvec':str(groupvec)})
            else:
                for p_nH, p_z, p_nMg, p_dG0 in pmap.ToMatrix():
                    rowdicts.append({'ID':mol_id, 'nH':p_nH, 'charge':p_z, 'nMg':p_nMg,
                                     'dG0':"%.1f" % p_dG0, 'groupvec':str(groupvec)})
        except GroupDecompositionError:
            rowdicts.append({'ID':mol_id, 'error':"cannot decompose"})
        except GroupMissingTrainDataError:
            rowdicts.append({'ID':mol_id, 'groupvec':str(groupvec),
                             'error':"missing training data"})
        
    if options.csv_output_filename is not None:
        out_fp = open(options.csv_output_filename, 'w')
        print "writing results to %s ... " % options.csv_output_filename
    else:
        out_fp = sys.stdout
    
    if options.biochemical:
        titles = ['ID', 'error', 'pH', 'I', 'pMg', 'dG0\'', 'groupvec']
    else:
        titles = ['ID', 'error', 'nH', 'nMg', 'charge', 'dG0', 'groupvec'] 
    csv_writer = csv.DictWriter(out_fp, titles)
    csv_writer.writeheader()
    csv_writer.writerows(rowdicts)
    
if __name__ == '__main__':
    CalculateThermo()