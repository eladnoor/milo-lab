import json, gzip, sys
from optparse import OptionParser
from pygibbs.kegg import Kegg
from pygibbs.nist_verify import LoadAllEstimators
from toolbox.database import SqliteDatabase
from pygibbs.kegg_compound import Compound
from pygibbs.kegg_reaction import Reaction
from pygibbs.kegg_enzyme import Enzyme

class KeggEncoder(json.JSONEncoder):
    def default(self, obj):
        if (isinstance(obj, Compound) or 
            isinstance(obj, Reaction) or 
            isinstance(obj, Enzyme)):
            return obj.ToJSONDict()
        return json.JSONEncoder.default(self, obj)

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
    opt_parser.add_option("-n", "--nullspace",
                          dest="nullspace_out_filename",
                          default="../res/kegg_gc_nullspace.json",
                          help="Null-space matrix for the GC estimations")
    opt_parser.add_option("-t", "--thermo_estimator",
                          dest="thermo_estimator",
                          default="PGC",
                          help="options are: " + ', '.join(estimators.keys()))
    return opt_parser

def WriteJSONFile(obj, filename):
    print 'Output filename: ', filename + '.gz'    
    fp = gzip.open(filename + '.gz', 'w')
    json.dump(obj, fp, cls=KeggEncoder, sort_keys=True, indent=4)
    fp.close()
    print 'Done'
    
def ExportJSONFiles():
    estimators = LoadAllEstimators()
    options, _ = MakeOpts(estimators).parse_args(sys.argv)
    
    print "Using the thermodynamic estimations of: " + options.thermo_estimator
    thermo = estimators[options.thermo_estimator]

    # Make sure we have all the data.
    kegg = Kegg.getInstance()
    kegg.AddThermodynamicData(estimators['alberty'], priority=1)
    kegg.AddThermodynamicData(thermo, priority=2)
    
    db = SqliteDatabase('../res/gibbs.sqlite')
    kegg.AddGroupVectorData(db, table_name='gc_groupvector')

    print 'Exporting Group Contribution Nullspace matrix as JSON.'
    nullspace_vectors = [row['group_vector'] 
                         for row in db.DictReader('gc_nullspace')]
    WriteJSONFile(nullspace_vectors, options.nullspace_out_filename)
        
    print 'Exporting KEGG compounds as JSON.'
    WriteJSONFile(kegg.AllCompounds(), options.compounds_out_filename)

    print 'Exporting KEGG reactions as JSON.'
    WriteJSONFile(kegg.AllReactions(), options.reactions_out_filename)
    
    print 'Exporting KEGG enzymes as JSON.'
    WriteJSONFile(kegg.AllEnzymes(), options.enzymes_out_filename)
    
if __name__ == '__main__':
    ExportJSONFiles()