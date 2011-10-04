#!/usr/bin/python

import json, gzip, sys
from optparse import OptionParser
from pygibbs.kegg import Kegg
from pygibbs.nist_verify import LoadAllEstimators
from toolbox.database import SqliteDatabase
from pygibbs.kegg_compound import Compound

class KeggEncoder(json.JSONEncoder):
    def default(self, obj):
        if (isinstance(obj, Compound)):
            return obj.ToJSONDict(light=True)
        return json.JSONEncoder.default(self, obj)

def MakeOpts(estimators):
    """Returns an OptionParser object with all the default options."""
    opt_parser = OptionParser()
    opt_parser.add_option("-o", "--out_filename",
                          dest="out_filename",
                          default="../res/kegg_pseudoisomers.json",
                          help="Output filename.")
    opt_parser.add_option("-s", "--thermodynamics_source",
                          dest="thermodynamics_source",
                          type="choice",
                          choices=estimators.keys(),
                          default="PGC",
                          help="The thermodynamic data to use")
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
    
    thermo = estimators[options.thermodynamics_source]
    print "Using the thermodynamic estimations of: " + thermo.name

    # Make sure we have all the data.
    kegg = Kegg.getInstance()
    kegg.AddThermodynamicData(estimators['alberty'], priority=1)
    kegg.AddThermodynamicData(thermo, priority=2)
    
    db = SqliteDatabase('../res/gibbs.sqlite')
    kegg.AddGroupVectorData(db, table_name='pgc_groupvector')

    print 'Exporting KEGG compound pseudoisomers as JSON.'
    WriteJSONFile(kegg.AllCompounds(), options.out_filename)

    
if __name__ == '__main__':
    ExportJSONFiles()