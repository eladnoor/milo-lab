#!/usr/bin/python

import logging
import sys
import sequence_utils

from Bio.Seq import Seq
from Bio.Seq import MutableSeq
from Bio.Alphabet import DNAAlphabet
from optparse import OptionParser


def MakeOpts():
    """Returns an OptionParser object with all the default options."""
    opt_parser = OptionParser()
    opt_parser.add_option("-l", "--left_sequence", dest="left_filename",
                          help="The filename of the left-side sequence.")
    opt_parser.add_option("-r", "--right_sequence",
                          dest="right_filename",
                          help="The filename of the left-side sequence.")
    opt_parser.add_option("-b", "--barcode_sequence",
                          dest="barcode_sequence",
                          help="The barcode to put between the two sequences")
    opt_parser.add_option("-o", "--overlap_length", type="int",
                          dest="overlap_length", default=20,
                          help="The amount of overlap between primers and each sequence.")
    
    return opt_parser


def main():
    """Makes primers for fusing two sequences using the Gibson method."""
    options, _ = MakeOpts().parse_args(sys.argv)
    
    if not options.left_filename or not options.right_filename:
        logging.fatal(
            'You must specify filenames for left and right sequences')
        
    print 'Left-side filename:', options.left_filename
    print 'Right-side filename:', options.right_filename
    
    alphabet = DNAAlphabet()
    left_sequence  = sequence_utils.MutableSeqFromFile(options.left_filename,
                                                       alphabet)
    right_sequence = sequence_utils.MutableSeqFromFile(options.right_filename,
                                                       alphabet)
    
    overlap = options.overlap_length
    primer = left_sequence[-overlap:]
    if options.barcode_sequence:
        bcode = Seq(options.barcode_sequence.upper(), alphabet)
        print 'Using barcode:', bcode
        primer += bcode
    primer += right_sequence[:overlap]
    
    print 'Primer length:', len(primer)
    print ('Forward primer (for right-side sequence "%s"):'
           % options.right_filename)
    print '\t%s' % primer
    primer.reverse_complement()
    print ('Reverse primer (for left-hand sequence "%s"):'
           % options.left_filename)
    print '\t%s' % primer
    

if __name__ == '__main__':
    main()