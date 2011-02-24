#!/usr/bin/python

import logging
import sys

from Bio.Seq import MutableSeq
from Bio.Seq import Seq
from Bio.Alphabet import DNAAlphabet
from optparse import OptionParser
import sequence_utils


ALPHABET = DNAAlphabet()
DEFAULT_SEQ_BASES = 20
START_HIS = Seq('ATGCATCATCACCATCACCAC', ALPHABET)
INV_CAP_LINKER = Seq('GCTAGCGTTGATCGGGCACGTAAGAG', ALPHABET)


def MakeOpts():
    """Returns an OptionParser object with all the default options."""
    opt_parser = OptionParser()
    opt_parser.add_option("-i", "--input_filename", dest="input_filename",
                          help="The filename of the sequence to make primers for.")
    opt_parser.add_option("-o", "--overlap_length", type="int",
                          dest="overlap_length", default=20,
                          help="The amount of overlap between the primer and sequence (BP).")    
    return opt_parser


def MakeForwardPrimer(seq, extension,
                      seq_bases=DEFAULT_SEQ_BASES):
    """Assumes the extension contains a start-codon."""
    my_seq = sequence_utils.DropStartCodon(seq)    
    return extension + my_seq[:seq_bases]


def MakeReversePrimer(seq, extension,
                      seq_bases=DEFAULT_SEQ_BASES):
    primer = seq[-seq_bases:] + extension
    primer.reverse_complement()
    return primer
    

def main():
    options, _ = MakeOpts().parse_args(sys.argv)

    if not options.input_filename:
        print 'Paste sequence you want to cut:',
        seq = MutableSeq(raw_input(), ALPHABET)
    else:
        seq = sequence_utils.MutableSeqFromFile(options.input_filename,
                                                ALPHABET)
    
    print 'Making primers with %d base-pairs sequence overlap' % options.overlap_length    
    fp = MakeForwardPrimer(seq, START_HIS, options.overlap_length)    
    print 'Forward primer (%d BP) with start codon & His-tag:' % len(fp)
    print fp
    

    rp = MakeReversePrimer(seq, INV_CAP_LINKER,
                           options.overlap_length)
    print 'Reverse primer (%d BP) with CAP linker:' % len(rp)
    print rp
    
    
if __name__ == '__main__':
    main()
    
    