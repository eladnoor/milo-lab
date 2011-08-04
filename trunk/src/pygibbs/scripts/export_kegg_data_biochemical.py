import csv, gzip, sys
from optparse import OptionParser
from pygibbs.kegg import Kegg
from pygibbs.nist_verify import LoadAllEstimators
from pygibbs.thermodynamic_errors import MissingCompoundFormationEnergy,\
    MissingReactionEnergy
from pygibbs.thermodynamic_constants import default_pH, default_I, default_T,\
    default_pMg
from pygibbs.kegg_errors import KeggReactionNotBalancedException,\
    KeggParseException

def MakeOpts(estimators):
    """Returns an OptionParser object with all the default options."""
    opt_parser = OptionParser()
    opt_parser.add_option("-H", "--ph", dest="pH", type='float', default=default_pH)
    opt_parser.add_option("-I", "--ionic_strength", type='float', dest="I", default=default_I)
    opt_parser.add_option("-M", "--pmg", dest="pMg", type='float', default=default_pMg)
    opt_parser.add_option("-T", "--temperature", dest="T", type='float', default=default_T)
    opt_parser.add_option("-c", "--compounds_out_filename",
                          dest="compounds_out_filename",
                          default="../res/kegg_compounds.csv",
                          help="Compounds output filename.")
    opt_parser.add_option("-r", "--reactions_out_filename",
                          dest="reactions_out_filename",
                          default="../res/kegg_reactions.csv",
                          help="Reactions output filename.")
    opt_parser.add_option("-t", "--thermo_estimator",
                          dest="thermo_estimator",
                          default="PGC",
                          help="options are: " + ', '.join(estimators.keys()))
    return opt_parser

def WriteCompoundCSV(compounds, thermo, filename):
    print 'Output filename: ', filename + '.gz'    
    fp = gzip.open(filename + '.gz', 'w')
    writer = csv.writer(fp)
    writer.writerow(["KEGG_CID", "dG0_tag", "pH", "I", "pMg", "T"])
    for compound in compounds:
        try:
            dG0_tag = "%.1f" % compound.PredictFormationEnergy(thermo)
        except MissingCompoundFormationEnergy as e:
            dG0_tag = None
        writer.writerow([compound.cid, dG0_tag, thermo.pH, thermo.I, 
                         thermo.pMg, thermo.T])
    fp.close()
    print 'Done'

def WriteReactionCSV(reactions, thermo, filename):
    print 'Output filename: ', filename + '.gz'    
    fp = gzip.open(filename + '.gz', 'w')
    writer = csv.writer(fp)
    writer.writerow(["KEGG_RID", "dG0_tag", "pH", "I", "pMg", "T"])
    for reaction in reactions:
        try:
            reaction.Balance(balance_water=True, exception_if_unknown=True)
            dG0_tag = "%.1f" % reaction.PredictReactionEnergy(thermo)
        except MissingCompoundFormationEnergy:
            dG0_tag = 'Missing formation energies'
        except KeggParseException as e:
            dG0_tag = str(e)
        except KeggReactionNotBalancedException:
            dG0_tag = 'Reaction cannot be balanced'
        except MissingReactionEnergy as e:
            dG0_tag = str(e)
        writer.writerow([reaction.rid, dG0_tag, thermo.pH, thermo.I,
                         thermo.pMg, thermo.T])
    fp.close()
    print 'Done'
    
def ExportCSVFiles():
    estimators = LoadAllEstimators()
    options, _ = MakeOpts(estimators).parse_args(sys.argv)
    
    print "Using the thermodynamic estimations of: " + options.thermo_estimator
    thermo = estimators[options.thermo_estimator]
    thermo.pH = float(options.pH)
    thermo.I = float(options.I)
    thermo.pMg = float(options.pMg)
    thermo.T = float(options.T)

    # Make sure we have all the data.
    kegg = Kegg.getInstance()
    
    print 'Exporting KEGG compounds as JSON.'
    WriteCompoundCSV(kegg.AllCompounds(), thermo, options.compounds_out_filename)

    print 'Exporting KEGG reactions as JSON.'
    WriteReactionCSV(kegg.AllReactions(), thermo, options.reactions_out_filename)
    
    
if __name__ == '__main__':
    ExportCSVFiles()