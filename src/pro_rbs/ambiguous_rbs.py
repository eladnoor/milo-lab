#!/usr/bin/python

import logging
import os
import pylab
import sys

from optparse import OptionParser
from pro_rbs import rbs_util
from scipy import stats
from toolbox import ambiguous_seq


# NOTE(flamholz): set this to the location of your NUPACK install.
os.environ['NUPACKHOME'] = '/home/flamholz/Downloads/nupack3.0/'


def MakeOpts():
    """Returns an OptionParser object with all the default options."""
    opt_parser = OptionParser()
    opt_parser.add_option("-r", "--rbs_seq", dest="rbs_seq",
                          default="../res/gibbs.sqlite",
                          help=("The RBS sequence including spacers, start codon, and"
                                "(optionally) his_tag. May be ambiguous."))
    opt_parser.add_option("-s", "--start_codon_seq", dest="start_codon_seq",
                          default="ATGCATCATCACCATCACCAC",
                          help="The sequence of the start codon, including some trailing bases.")
    return opt_parser


def HandleAmbigSeq(ambig_seq, start_codon_seq):
    """Process an ambiguous RBS sequence.
    
    Test all concrete possibilities for the sequence, draw a histogram of 
    their binding energies.
    
    Args:
        ambig_seq: a AmbigousDNASeq object.
        start_codon_seq: the sequence of the start codon including
            some trailing bases.
    """
    dg_by_seq = {}
    failure_count = 0
    for concrete_seq in ambig_seq.AllConcreteSeqs():
        try:
            seq_str = concrete_seq.tostring()
            start_index = seq_str.find(start_codon_seq)
            if start_index == -1:
                logging.error('Couldn\'t find start codon.')
                continue
            
            dG = rbs_util.TryCalcDG(concrete_seq,
                                    start_codon_index=start_index)
            dg_by_seq[str(concrete_seq)] = dG
        except Exception, e:
            logging.error(e)
            failure_count += 1
            continue
    
    print 'Tested %d RBS options.' % len(dg_by_seq)
    print '%d RBS options failed.' % failure_count
    
    binding_dGs = dg_by_seq.values()
    log_rates = map(rbs_util.LogRate, binding_dGs)
    kde = stats.gaussian_kde(log_rates)
    sample_points = pylab.linspace(-4, 4, 1000)
    kde_pts = kde.evaluate(sample_points)

    pylab.xlim((-4, 4))    
    pylab.xlabel('log10(K_binding)')
    pylab.ylabel('Number of RBS')
    pylab.hist(log_rates, label='RBS histogram')
    pylab.plot(sample_points, kde_pts, color='r', label='RBS KDE')
    pylab.legend()
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
    start_codon_seq = options.start_codon_seq.upper()
    print 'Start codon sequence:', start_codon_seq
    
    ambig_seq = ambiguous_seq.AmbigousDNASeq(options.rbs_seq)
    if ambig_seq.IsConcrete():
        print 'Got a concrete sequence. Bailing.'
        return
    
    print 'Testing all options for ambiguous sequence'
    HandleAmbigSeq(ambig_seq, start_codon_seq)
        

if __name__ == '__main__':
    main()
    