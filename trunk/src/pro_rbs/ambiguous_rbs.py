#!/usr/bin/python

import logging
import os
import pylab
import sys

from optparse import OptionParser
from pro_rbs import rbs_util
from toolbox import ambiguous_seq


# NOTE(flamholz): set this to the location of your NUPACK install.
os.environ['NUPACKHOME'] = '/home/flamholz/Downloads/nupack3.0/'


def MakeOpts():
    """Returns an OptionParser object with all the default options."""
    opt_parser = OptionParser()
    opt_parser.add_option("-s", "--rbs_seq", dest="rbs_seq",
                          default="../res/gibbs.sqlite",
                          help=("The RBS sequence including spacers, start codon, and"
                                "(optionally) his_tag. May be ambiguous."))
    return opt_parser


def HandleAmbigSeq(ambig_seq):
    """Process an ambiguous RBS sequence.
    
    Test all concrete possibilities for the sequence, draw a histogram of 
    their binding energies.
    
    Args:
        ambig_seq: a AmbigousDNASeq object.
    """
    dg_by_seq = {}
    failure_count = 0
    for concrete_seq in ambig_seq.AllConcreteSeqs():
        try:
            dG = rbs_util.TryCalcDG(concrete_seq)
            dg_by_seq[str(concrete_seq)] = dG
        except Exception, e:
            logging.error(e)
            failure_count += 1
            continue
    
    print 'Tested %d RBS options.' % len(dg_by_seq)
    print '%d RBS options failed.' % failure_count
    
    binding_dGs = dg_by_seq.values()
    log_rates = map(rbs_util.LogRate, binding_dGs)
    pylab.hist(log_rates)
    pylab.xlim((-4, 4))
    pylab.xlabel('log10(K_binding)')
    pylab.ylabel('Number of RBS')
    pylab.show()
    
    seq_by_dg = {}
    for seq, dg in dg_by_seq.iteritems():
        rounded = round(dg, 1)
        seq_by_dg.setdefault(rounded, []).append(seq)
        
    print 'sequence by dG'

    for dG in sorted(seq_by_dg.keys()):
        print '%.2g:' % dG,
        for seq in seq_by_dg[dG]:
            print seq


def main():
    options, _ = MakeOpts().parse_args(sys.argv)
    print 'RBS sequence:', options.rbs_seq
    
    ambig_seq = ambiguous_seq.AmbigousDNASeq(options.rbs_seq)
    if ambig_seq.IsConcrete():
        print 'Got a concrete sequence. Bailing.'
        return
    
    print 'Testing all options for ambiguous sequence'
    HandleAmbigSeq(ambig_seq)
        

if __name__ == '__main__':
    main()
    