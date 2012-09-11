import sys
from optparse import OptionParser
from pygibbs.nist_verify import LoadAllEstimators
from toolbox.database import SqliteDatabase
from pygibbs.groups import GroupContribution
from toolbox.html_writer import NullHtmlWriter
from pygibbs.thermodynamic_errors import MissingCompoundFormationEnergy,\
    MissingReactionEnergy
import csv
from pygibbs.kegg_errors import KeggReactionNotBalancedException,\
    KeggParseException
from pygibbs.thermodynamic_constants import default_pH, default_pMg, default_I,\
    default_T
from toolbox.molecule import OpenBabelError

def MakeOpts(estimators):
    """Returns an OptionParser object with all the default args."""
    opt_parser = OptionParser()
    opt_parser.add_option("-c", "--compounds_out_filename",
                          dest="compounds_out_filename",
                          type="string", action="store",
                          default="../res/pgc_compounds.csv",
                          help="Compounds output filename")
    opt_parser.add_option("-r", "--reactions_out_filename",
                          dest="reactions_out_filename",
                          default="../res/pgc_reactions.csv",
                          help="Reactions output filename")
    opt_parser.add_option("-b", "--biochemical",
                          dest="transformed",
                          action="store_true", default=False,
                          help="Use biochemical (transformed) Group Contributions")
    opt_parser.add_option("-p", "--ph", dest="ph", type="float",
                          default=default_pH, help="The pH")
    opt_parser.add_option("-m", "--pmg", dest="pmg", type="float",
                          default=default_pMg, help="The pMg")
    opt_parser.add_option("-i", "--ionic_strength", dest="i_s", type="float",
                          default=default_I, help="The ionic strength, in M")
    opt_parser.add_option("-t", "--temperature", dest="temp", type="float",
                          default=default_T, help="The temperature, in K")
    return opt_parser

def main():
    estimators = LoadAllEstimators()
    args, _ = MakeOpts(estimators).parse_args(sys.argv)
    
    # Make sure we have all the data.
    db = SqliteDatabase('../res/gibbs.sqlite')
    G = GroupContribution(db=db, html_writer=NullHtmlWriter(),
                          transformed=args.transformed)
    G.init()
    
    print 'Exporting KEGG compounds to %s' % args.compounds_out_filename
    csv_writer = csv.writer(open(args.compounds_out_filename, 'w'))
    csv_writer.writerow(["KEGG ID", "nH", "CHARGE", "nMg", "dG0_f"])
    for cid in sorted(G.get_all_cids()):
        try:
            for nH, z, nMg, dG0 in G.cid2PseudoisomerMap(cid).ToMatrix():
                csv_writer.writerow(["C%05d" % cid, nH, z, nMg, "%.1f" % dG0])
        except MissingCompoundFormationEnergy as e:
            csv_writer.writerow(["C%05d" % cid, None, None, None, str(e)])
        
    print 'Exporting KEGG reactions to %s' % args.reactions_out_filename
    csv_writer = csv.writer(open(args.reactions_out_filename, 'w'))
    csv_writer.writerow(["KEGG ID", "dG'0_r (pH=%.1f, I=%.2f, pMg=%.1f, T=%.1f)" % 
                         (args.ph, args.i_s, args.pmg, args.temp)])
    for rid in sorted(G.kegg.get_all_rids()):
        reaction = G.kegg.rid2reaction(rid)
        try:
            reaction.Balance(balance_water=True)
            dG0_r = reaction.PredictReactionEnergy(G, pH=args.ph,
                        pMg=args.pmg, I=args.i_s, T=args.temp)
            csv_writer.writerow(["R%05d" % rid, "%.1f" % dG0_r])
        except (KeggParseException,
                MissingCompoundFormationEnergy, 
                KeggReactionNotBalancedException,
                MissingReactionEnergy,
                KeyError,
                OpenBabelError) as e:
            csv_writer.writerow(["R%05d" % rid, str(e)])
    
if __name__ == '__main__':
    main()