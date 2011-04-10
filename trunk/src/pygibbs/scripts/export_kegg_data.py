
import json
import sys

from optparse import OptionParser
from pygibbs.kegg import Kegg
from toolbox.database import SqliteDatabase
from pygibbs.thermodynamics import PsuedoisomerTableThermodynamics


def MakeOpts():
    """Returns an OptionParser object with all the default options."""
    opt_parser = OptionParser()
    opt_parser.add_option("-k", "--kegg_database_location", 
                          dest="kegg_db_filename",
                          default="../data/public_data.sqlite",
                          help="The KEGG database location")
    opt_parser.add_option("-c", "--compounds_out_filename",
                          dest="compounds_out_filename",
                          default="../res/kegg_compounds.json",
                          help="Compounds output filename.")
    opt_parser.add_option("-r", "--reactions_out_filename",
                          dest="reactions_out_filename",
                          default="../res/kegg_reactions.json",
                          help="Reactions output filename.")
    opt_parser.add_option("-e", "--enzymes_out_filename",
                          dest="enzymes_out_filename",
                          default="../res/kegg_enzymes.json",
                          help="Enzymes output filename.")
    opt_parser.add_option("-t", "--thermodynamics_filename",
                          dest="thermo_filename",
                          default='../data/thermodynamics/dG0.csv',
                          help="Thermodynamics filename")
    opt_parser.add_option("-d", "--database_location", 
                          dest="db_filename",
                          default="../res/gibbs.sqlite",
                          help="The Thermodynamic database location")
    opt_parser.add_option("-g", "--gc_table_name",
                          dest="gc_table_name",
                          default='gc_pseudoisomers',
                          help="Group Contribution Table Name")    
    return opt_parser


def WriteJSONFile(objects, filename):
    json_dicts = [item.ToJSONDict() for item in objects]
    json_file = open(filename, 'w')
    json.dump(json_dicts, json_file, sort_keys=True, indent=4)
    json_file.close()
    

def ExportJSONFiles():
    options, _ = MakeOpts().parse_args(sys.argv)
    print 'KEGG Database filename:', options.kegg_db_filename
    print 'Observed Thermodynamics filename:', options.thermo_filename
    print 'Thermodynamic Database filename:', options.db_filename
    print 'Group Contribution Table Name:', options.gc_table_name

    db = SqliteDatabase(options.db_filename)
    observed_thermo = PsuedoisomerTableThermodynamics.FromCsvFile(
        options.thermo_filename)
    if not db.DoesTableExist(options.gc_table_name):
        raise ValueError('The table %s does not exist in the database. '
                         'Please run the groups.py script and try again.'
                         % options.gc_table_name)
    thermo = PsuedoisomerTableThermodynamics.FromDatabase(
        db, options.gc_table_name)
    thermo.override_data(observed_thermo)

    # Make sure we have all the data.
    kegg = Kegg.getInstance()
    kegg.AddThermodynamicData(thermo)
    
    print 'Exporting KEGG compounds as JSON.'
    print 'Output filename:', options.compounds_out_filename    
    WriteJSONFile(kegg.AllCompounds(),
                  options.compounds_out_filename)
    print 'Done'

    print 'Exporting KEGG reactions as JSON.'
    print 'Output filename:', options.reactions_out_filename
    WriteJSONFile(kegg.AllReactions(),
                  options.reactions_out_filename)
    print 'Done'
    
    print 'Exporting KEGG enzymes as JSON.'
    print 'Output filename:', options.enzymes_out_filename
    WriteJSONFile(kegg.AllEnzymes(),
                  options.enzymes_out_filename)
    print 'Done'
    

if __name__ == '__main__':
    ExportJSONFiles()