import json
import sys
from optparse import OptionParser
from pygibbs.kegg import Kegg
from pygibbs.nist_verify import LoadAllEstimators

def MakeOpts(estimators):
    """Returns an OptionParser object with all the default options."""
    opt_parser = OptionParser()
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
    opt_parser.add_option("-t", "--thermo_estimator",
                          dest="thermo_estimator",
                          default="milo_gc",
                          help="options are: " + ', '.join(estimators.keys()))
    return opt_parser


def WriteJSONFile(objects, filename):
    json_dicts = [item.ToJSONDict() for item in objects]
    json_file = open(filename, 'w')
    json.dump(json_dicts, json_file, sort_keys=True, indent=4)
    json_file.close()
    

def ExportJSONFiles():
    estimators = LoadAllEstimators()
    options, _ = MakeOpts(estimators).parse_args(sys.argv)
    
    print "Using the thermodynamic estimations of: " + options.thermo_estimator
    thermo = estimators[options.thermo_estimator]

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