#!/usr/bin/python

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

# Headers
KEGG_REACTION = '!MiriamID::urn:miriam:kegg.reaction'
KEGG_COMPOUND = '!MiriamID::urn:miriam:kegg.compound'
NAME = '!Name'
FORMATION_ENERGY = '!dG0_tag (kJ/mol)'
REACTION_ENERGY = '!dG0_tag (kJ/mol)'
COMPOUND_ROW_ORDER = [KEGG_COMPOUND, NAME, FORMATION_ENERGY,
                      '!pH', '!I (mM)', '!T (Kelvin)', '!Note']
REACTION_ROW_ORDER = [KEGG_REACTION, REACTION_ENERGY,
                      '!pH', '!I (mM)', '!T (Kelvin)', '!Note']


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
    fname = filename + '.gz'
    print 'Output filename: ', fname    
    fp = gzip.open(fname, 'w')
    writer = csv.writer(fp)
    writer.writerow(COMPOUND_ROW_ORDER)
    for compound in compounds:
        dG0_tag, note = None, None
        try:
            dG0_tag = "%.1f" % compound.PredictFormationEnergy(thermo)
        except MissingCompoundFormationEnergy, e:
            note = str(e)
        writer.writerow([compound.cid, compound.name, dG0_tag, thermo.pH, thermo.I, 
                         thermo.T, note])
    fp.close()
    print 'Done'


def WriteReactionCSV(reactions, thermo, filename):
    fname = filename + '.gz'
    print 'Output filename: ', fname    
    fp = gzip.open(fname, 'w')
    writer = csv.writer(fp)
    writer.writerow(REACTION_ROW_ORDER)
    for reaction in reactions:
        dG0_tag, note = None, None
        try:
            reaction.Balance(balance_water=True, exception_if_unknown=True)
            dG0_tag = "%.1f" % reaction.PredictReactionEnergy(thermo)
        except MissingCompoundFormationEnergy:
            note = 'Missing formation energies'
        except KeggParseException as e:
            note = str(e)
        except KeggReactionNotBalancedException, e:
            note = str(e)
        except MissingReactionEnergy as e:
            note = str(e)
        writer.writerow([reaction.rid, dG0_tag, thermo.pH, thermo.I,
                         thermo.T, note])
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