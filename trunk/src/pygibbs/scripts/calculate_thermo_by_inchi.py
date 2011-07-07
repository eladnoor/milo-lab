import sys
import csv
from optparse import OptionParser
from pygibbs.groups import GroupContribution, GroupMissingTrainDataError
from toolbox.database import SqliteDatabase
from toolbox.molecule import Molecule
from pygibbs.group_decomposition import GroupDecompositionError

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
        inchi = row["InChI"]
        try:
            mol = Molecule.FromInChI(inchi)
            decomposition = G.Mol2Decomposition(mol, 
                ignore_protonations=ignore_protonations)
            nH = decomposition.Hydrogens()
            z = decomposition.NetCharge()
            nMg = decomposition.Magnesiums()
            groupvec = decomposition.AsVector()
            dG0 = G.groupvec2val(groupvec)
            if nMg > 0:
                raise Exception("Magnesium ions are not allowed here")
            rowdicts.append({'ID':id, 'nH':nH, 'charge':z, 'dG0':dG0,
                              'groupvec':str(groupvec), 'error':None})
        except GroupDecompositionError:
            rowdicts.append({'ID':id, 'nH':None, 'charge':None, 'dG0':None,
                              'groupvec':None, 'error':"cannot decompose"})
        except GroupMissingTrainDataError:
            rowdicts.append({'ID':id, 'nH':nH, 'charge':z, 'dG0':None, 
                              'groupvec':str(groupvec), 'error':"missing training data"})
        
    print "writing results to %s" % options.csv_output_filename
    csv_writer = csv.DictWriter(open(options.csv_output_filename, 'w'), 
                                ['ID', 'error', 'nH', 'charge', 'dG0', 'groupvec'])
    csv_writer.writeheader()
    csv_writer.writerows(rowdicts)
    
if __name__ == '__main__':
    CalculateThermo()