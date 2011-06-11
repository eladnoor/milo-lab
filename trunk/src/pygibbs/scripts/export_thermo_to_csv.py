from pygibbs.thermodynamics import PsuedoisomerTableThermodynamics
from pygibbs.pseudoisomer import PseudoisomerMap
from optparse import OptionParser
from toolbox.database import SqliteDatabase
import sys

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
    opt_parser.add_option("-g", "--thermo_table_name",
                          dest="thermo_table_name",
                          default='hatzi_thermodynamics',
                          help="Thermodynamics Table Name")    
    opt_parser.add_option("-c", "--csv_out_filename",
                          dest="csv_out_filename",
                          default="../res/thermo.csv",
                          help="CSV output filename")
    return opt_parser

def ExportThermo():
    options, _ = MakeOpts().parse_args(sys.argv)
    print 'KEGG Database filename:', options.kegg_db_filename
    print 'Thermodynamic Database filename:', options.db_filename
    print 'Group Contribution Table Name:', options.thermo_table_name
    print 'CSV output filename:', options.csv_out_filename

    db = SqliteDatabase(options.db_filename)
    thermo = PsuedoisomerTableThermodynamics.FromDatabase(db, 'hatzi_thermodynamics')
    thermo.AddPseudoisomer( 139, nH=0,  z=1, nMg=0, dG0=0)      # Ferrodoxin(ox)
    thermo.AddPseudoisomer( 138, nH=0,  z=0, nMg=0, dG0=38.0)   # Ferrodoxin(red)
    thermo.AddPseudoisomer( 399, nH=90, z=0, nMg=0, dG0=0)      # Ubiquinone-10(ox)
    thermo.AddPseudoisomer( 390, nH=92, z=0, nMg=0, dG0=-103.2) # Ubiquinone-10(red)
    thermo.AddPseudoisomer( 828, nH=16, z=0, nMg=0, dG0=0)      # Menaquinone(ox)
    thermo.AddPseudoisomer(5819, nH=18, z=0, nMg=0, dG0=-65.8)  # Menaquinone(red)
    thermo.SetPseudoisomerMap(101, PseudoisomerMap(nH=23, z=0, nMg=0, dG0=0.0)) # THF
    thermo.SetPseudoisomerMap(234, PseudoisomerMap(nH=23, z=0, nMg=0, dG0=-137.5)) # 10-Formyl-THF
    thermo.SetPseudoisomerMap(445, PseudoisomerMap(nH=22, z=0, nMg=0, dG0=65.1)) # 5,10-Methenyl-THF
    thermo.SetPseudoisomerMap(143, PseudoisomerMap(nH=23, z=0, nMg=0, dG0=77.9)) # 5,10-Methylene-THF
    thermo.SetPseudoisomerMap(440, PseudoisomerMap(nH=25, z=0, nMg=0, dG0=32.1)) # 5-Methyl-THF
    
    thermo.write_data_to_csv(options.csv_out_filename)
    
if __name__ == "__main__":
    ExportThermo()
