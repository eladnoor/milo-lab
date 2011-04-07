#!/usr/bin/python

import logging
import numpy
import os
import pylab
import sys

from optparse import OptionParser
from pro_rbs import rbs_util
from scipy import stats
from toolbox import ambiguous_seq
from toolbox import random_seq


# NOTE(flamholz): set this to the location of your NUPACK install.
os.environ['NUPACKHOME'] = '/home/flamholz/Downloads/nupack3.0/'


def MakeOpts():
    """Returns an OptionParser object with all the default options."""
    opt_parser = OptionParser()
    opt_parser.add_option("-r", "--rbs_seq", dest="rbs_seq",
                          help=("The RBS sequence including spacers, start codon, and"
                                "(optionally) his_tag. Cannot be ambiguous."))
    opt_parser.add_option("-s", "--start_codon_seq", dest="start_codon_seq",
                          default="ATGCATCATCACCATCACCAC",
                          help="The sequence of the start codon, including some trailing bases.")
    opt_parser.add_option("-n", "--num_trials", dest="num_trials", type="int",
                          default=100,
                          help="The number of trials")
    return opt_parser


def HandleSeq(seq, start_codon_seq, num_trials):
    """Process an RBS sequence.
    
    Test many random genes and see how much the RBS binding changes.
    
    Args:
        seq: a DNA sequence for the RBS.
        start_codon_seq: the sequence of the start codon including
            some trailing bases.
        num_trials: the number of random genes to test.
    """
    start_index = seq.find(start_codon_seq)
    if start_index == -1:
        logging.error('Couldn\'t find start codon.')
        return
    reference_dG = rbs_util.TryCalcDG(seq,
                                      start_codon_index=start_index)
    reference_log_rate = rbs_util.LogRate(reference_dG)
    
    dGs = []
    failure_count = 0
    for i in xrange(num_trials):
        gene_str = random_seq.RandomSeq(100).tostring()
        if i % 25 == 0:
            print 'Testing random gene sequence', i, gene_str
    
        try:
            full_rbs_str = seq + gene_str
            dG = rbs_util.TryCalcDG(full_rbs_str,
                                    start_codon_index=start_index)
            
            dGs.append(dG)
        except Exception, e:
            logging.error(e)
            failure_count += 1
            continue
        
    print 'Tested %d RBS options.' % len(dGs)
    print '%d RBS options failed.' % failure_count
    
    log_rates = map(rbs_util.LogRate, dGs)
    centered_log_rates = [rate - reference_log_rate for rate in log_rates]
    print 'min(log(k_binding)) =', min(log_rates)
    print 'max(log(k_binding)) =', max(log_rates)
    print 'std(log(k_binding)) =', pylab.std(log_rates)
    
    pylab.xlim((-4, 4))
    pylab.xlabel('Reference minus log10(K_binding)')
    pylab.ylabel('Number of Sequences')
    pylab.hist(centered_log_rates, label='Binding Histogram for %s trials' % num_trials)
    pylab.legend()
    pylab.show()


def main():
    options, _ = MakeOpts().parse_args(sys.argv)
    print 'RBS sequence:', options.rbs_seq
    start_codon_seq = options.start_codon_seq.upper()
    print 'Start codon sequence:', start_codon_seq
    
    num_trials = options.num_trials
    print 'Running', num_trials, 'trials'
    
    HandleSeq(options.rbs_seq, start_codon_seq, num_trials)
        

if __name__ == '__main__':
    main()
    