#!/usr/bin/python

import sys
import pylab
from pygibbs.thermodynamic_constants import default_T, default_I, default_pMg
from pygibbs.kegg import Kegg
from pygibbs.nist import NistRowData, Nist
from pygibbs.kegg_errors import KeggParseException
from optparse import OptionParser
from pygibbs.nist_verify import LoadAllEstimators

def MakeOpts():
    """Returns an OptionParser object with all the default options."""
    usage = "usage: %prog [options] reaction"
    opt_parser = OptionParser(usage=usage)
    opt_parser.add_option("-r", "--reaction", action="store", type="int",
                          dest="rid", default=None,
                          help="The KEGG ID of a reaction")
    opt_parser.add_option("-p", "--ph", action="store", type="float",
                          dest="pH", default=2,
                          help="The width of the range of pH values to plot")
    opt_parser.add_option("-i", "--ionic_strength", action="store", type="float",
                          dest="I", default=default_I,
                          help="The Ionic strength of the solution (in M)")
    opt_parser.add_option("-m", "--pmg", action="store", type="float",
                          dest="pMg", default=default_pMg,
                          help="The pMg of the solution")
    opt_parser.add_option("-t", "--temperature", action="store", type="float",
                          dest="T", default=default_T,
                          help="The T of the solution (in K)")
    opt_parser.add_option("-o", "--output", action="store", type="string",
                          dest="output", default=None,
                          help="Destination filename for writing the figure")
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


def main():
    kegg = Kegg.getInstance()
    options, args = MakeOpts().parse_args(sys.argv)
    print ('Parameters: T=%f K, pMg=%.2g, I=%.2gM' % 
           (options.T, options.pMg, options.I))
    print "reaction:", args[-1]

    estimators = LoadAllEstimators()
    del estimators['alberty']
    del estimators['nist_regression']
    del estimators['hatzi_gc_pka']
    
    pylab.rcParams['text.usetex'] = False
    pylab.rcParams['legend.fontsize'] = 8
    pylab.rcParams['font.family'] = 'sans-serif'
    pylab.rcParams['font.size'] = 12
    pylab.rcParams['lines.linewidth'] = 2
    pylab.rcParams['lines.markersize'] = 10
    pylab.rcParams['lines.markeredgewidth'] = 0
    pylab.rcParams['figure.figsize'] = [5.0, 5.0]
    pylab.rcParams['figure.dpi'] = 100
    
    colors = {}
    colors['nist'] = [64.0/255, 111.0/255, 29.0/255, 1.0] # greenish
    colors['hatzi_gc'] = [54.0/255, 182.0/255, 202.0/255, 1.0] # cyanish
    colors['milo_gc'] = [202.0/255, 101.0/255, 54.0/255, 1.0] # orangish
    
    fig = pylab.figure()
    
    fig.hold(True)
    reaction = GetSparseReactionInput(args[-1], kegg)
    reaction.Balance()
    print 'Reaction: %s' % reaction.FullReactionString()

    nist = Nist()
    nist_rows = nist.SelectRowsFromNist(reaction, check_reverse=True)

    pH_min = 7.0 - options.pH/2.0
    pH_max = 7.0 + options.pH/2.0
    
    if nist_rows:
        dG0_list = []
        pH_list = []
        for row_data in nist_rows:
            pH_list.append(row_data.pH)
            if row_data.reaction == reaction:
                dG0_list.append(row_data.dG0_r)
            else:
                dG0_list.append(-row_data.dG0_r)
    
        pylab.plot(pH_list, dG0_list, marker='.', linestyle='None',
                   label='measured data', markerfacecolor=colors['nist'])
        pH_max = max(pH_list + [pH_max])
        pH_min = min(pH_list + [pH_min])
    
    pH_range = pylab.arange(pH_min-0.1, pH_max+0.1, 0.02)
    for key, thermo in estimators.iteritems():
        print key, 'dG0 at pH=7: %.2f' % reaction.PredictReactionEnergy(thermo, 
                pH=7.0, pMg=options.pMg, I=options.I, T=options.T)
        dG0 = []
        for pH in pH_range:
            dG0.append(reaction.PredictReactionEnergy(thermo, 
                pH=pH, pMg=options.pMg, I=options.I, T=options.T))
        pylab.plot(pH_range, dG0, marker='None', linestyle='solid', color=colors[key],
                   figure=fig, label=thermo.name)

    pylab.xlabel('pH')
    pylab.ylabel(r'$\Delta_r G^\circ$')
    pylab.title(kegg.reaction2string(reaction), fontsize=8)
    pylab.legend(loc='lower left')

    if not options.output:
        pylab.show()
    else:
        fig.savefig(options.output, format='svg')

if __name__ == '__main__':
    main()
    
