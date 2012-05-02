import sys
from optparse import OptionParser
from pygibbs.thermodynamic_constants import default_I, default_pH, default_pMg,\
    default_T
from pygibbs.nist_verify import LoadAllEstimators
from pygibbs.kegg_errors import KeggParseException
from pygibbs.nist import NistRowData
from pygibbs.kegg import Kegg

def MakeOpts(estimators):
    """Returns an OptionParser object with all the default options."""
    usage = "usage: %prog [options] reaction"
    opt_parser = OptionParser(usage=usage)
    opt_parser.add_option("-r", "--reaction", action="store", type="int",
                          dest="rid", default=None,
                          help="The KEGG ID of a reaction")
    opt_parser.add_option("-p", "--ph", action="store", type="float",
                          dest="pH", default=default_pH,
                          help="The pH of the solution")
    opt_parser.add_option("-i", "--ionic_strength", action="store", type="float",
                          dest="I", default=default_I,
                          help="The Ionic strength of the solution (in M)")
    opt_parser.add_option("-m", "--pmg", action="store", type="float",
                          dest="pMg", default=default_pMg,
                          help="The pMg of the solution")
    opt_parser.add_option("-t", "--temperature", action="store", type="float",
                          dest="T", default=default_T,
                          help="The Temperature of the solution (in K)")
    opt_parser.add_option("-s", "--thermodynamics_source",
                          dest="thermodynamics_source",
                          type="choice",
                          choices=estimators.keys(),
                          default="UGC",
                          help="The thermodynamic data to use")
    return opt_parser

def GetSparseReactionInput(args, kegg):
    if not args:
        sys.exit(0)
    elif args[0] == 'R':
        return kegg.rid2reaction(int(args[1:]))
    else:
        try:
            return NistRowData.ParseReactionFormula("a reaction", args)
        except KeggParseException:
            raise ValueError("Invalid reaction: " + args)

def CalculateThermo():
    estimators = LoadAllEstimators()
    parser = MakeOpts(estimators)
    options, args = parser.parse_args(sys.argv)

    kegg = Kegg.getInstance()
    if options.rid is None:
        reaction = GetSparseReactionInput(args[-1], kegg)
    else:
        reaction = kegg.rid2reaction(options.rid)
    
    estimator = estimators[options.thermodynamics_source]
    pH, I, pMg, T = options.pH, options.I, options.pMg, options.T
    estimator.SetConditions(pH=pH, I=I, pMg=pMg, T=T)

    print "Thermodynamic source:", options.thermodynamics_source
    print ('Parameters: pH=%.1f, pMg=%.1f, I=%.2fM, T=%.1fK' % 
           (options.pH, options.pMg, options.I, options.T))
    print str(reaction)
    print 'dG\'0 = %.2f [kJ/mol]' % reaction.PredictReactionEnergy(estimator)

if __name__ == '__main__':
    CalculateThermo()