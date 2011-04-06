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
                                "(optionally) his_tag. May be ambiguous."))
    opt_parser.add_option("-a", "--gene_a", dest="gene_a",
                          help="Sequence of a first gene.")
    opt_parser.add_option("-b", "--gene_b", dest="gene_b",
                          help="Sequence of a second gene.")
    opt_parser.add_option("-s", "--start_codon_seq", dest="start_codon_seq",
                          default="ATGCATCATCACCATCACCAC",
                          help="The sequence of the start codon, including some trailing bases.")
    return opt_parser


def HandleAmbigSeq(ambig_seq, start_codon_seq, gene_a, gene_b):
    """Process an ambiguous RBS sequence.
    
    Test all concrete possibilities for the sequence, draw a histogram of 
    their binding energies.
    
    Args:
        ambig_seq: a AmbigousDNASeq object for the RBS.
        start_codon_seq: the sequence of the start codon including
            some trailing bases.
        gene_a: the sequence of the first gene.
        gene_b: the sequence of the second gene.
    """
    dGs_a = []
    dGs_b = []
    failure_count = 0
    for i, concrete_seq in enumerate(ambig_seq.AllConcreteSeqs()):
        if i % 25 == 0:
            print 'Testing concrete sequence', i, concrete_seq
        
        try:
            seq_str = concrete_seq.tostring()
            start_index = seq_str.find(start_codon_seq)
            if start_index == -1:
                logging.error('Couldn\'t find start codon.')
                continue
            
            seq_a = seq_str + gene_a
            seq_b = seq_str + gene_b
            
            dG_a = rbs_util.TryCalcDG(seq_a,
                                      start_codon_index=start_index)
            dG_b = rbs_util.TryCalcDG(seq_b,
                                      start_codon_index=start_index)
            
            dGs_a.append(dG_a)
            dGs_b.append(dG_b)
        except Exception, e:
            logging.error(e)
            failure_count += 1
            continue
    
    
    assert len(dGs_a) == len(dGs_b)
    print 'Tested %d RBS options.' % len(dGs_a)
    print '%d RBS options failed.' % failure_count
    
    log_rates_a = map(rbs_util.LogRate, dGs_a)
    log_rates_b = map(rbs_util.LogRate, dGs_b)
    
    corr = numpy.corrcoef(log_rates_a, log_rates_b, rowvar=True)
    r = corr[0][1]
    print 'correlations coeffs', corr
    print 'R =', r
    
    pylab.xlabel('log10(K_binding) with gene A')
    pylab.ylabel('log10(K_binding) with gene A')
    pylab.scatter(log_rates_a, log_rates_b, label='gene A-gene B (r=%.2g)' % r)
    
    pylab.legend()
    pylab.show()


def main():
    options, _ = MakeOpts().parse_args(sys.argv)
    print 'RBS sequence:', options.rbs_seq
    start_codon_seq = options.start_codon_seq.upper()
    print 'Start codon sequence:', start_codon_seq
    
    gene_a = options.gene_a
    gene_b = options.gene_b
    
    if gene_a:
        gene_a = gene_a.upper()
    else:
        print 'No Gene A provided. Generating randomized gene.'
        gene_a = random_seq.RandomSeq(100).tostring().upper()

    if gene_b:
        gene_b = gene_b.upper()
    else:
        print 'No Gene B provided. Generating randomized gene.'
        gene_b = random_seq.RandomSeq(100).tostring()
        
    print 'Gene A:', gene_a
    print 'Gene B:', gene_b 
    
    ambig_seq = ambiguous_seq.AmbigousDNASeq(options.rbs_seq)
    if ambig_seq.IsConcrete():
        print 'Got a concrete sequence. Bailing.'
        return
    
    print 'Testing all options for ambiguous sequence'
    HandleAmbigSeq(ambig_seq, start_codon_seq, gene_a, gene_b)
        

if __name__ == '__main__':
    main()
    