#!/usr/bin/python

import getopt
import sys


DEFAULT_SEQ_BASES = 20
START_HIS = 'ATGCATCATCACCATCACCAC'
INV_CAP_LINKER = 'GCTAGCGTTGATCGGGCACGTAAGAG'

INVERSES = {'A': 'T',
            'T': 'A',
            'G': 'C',
            'C': 'G'}


def InvertAndReverse(seq):
    f = lambda b: INVERSES[b]
    l = map(f, seq)
    l.reverse()
    return ''.join(l)


def DropStartCodon(seq):
    if seq.startswith('ATG'):
        return seq[3:]
    return seq


def MakeForwardPrimer(seq, extension,
                      drop_start_codon=True,
                      seq_bases=DEFAULT_SEQ_BASES):
    if drop_start_codon:
        seq = DropStartCodon(seq)
        
    return '%s%s' % (extension, seq[:seq_bases])


def MakeReversePrimer(seq, extension,
                      seq_bases=DEFAULT_SEQ_BASES):
    primer = seq[-seq_bases:] + extension
    return InvertAndReverse(primer)
    

def PrintPrimers(seq, forward_extension, reverse_extension,
                 seq_bases=DEFAULT_SEQ_BASES, drop_start_codon=True):
    s = seq.upper()
    print 'Forward primer:', MakeForwardPrimer(s, forward_extension,
                                               drop_start_codon, seq_bases)
    print 'Reverse primer:', MakeReversePrimer(s, reverse_extension,
                                               seq_bases)
    

def main(argv=None):
    if argv is None:
        argv = sys.argv
        
    try:
        opts, unused_args = getopt.getopt(argv[1:], 'i:',
                                          ['input='])
    except getopt.error:
        print 'Failed to parse opt'
        sys.exit()
    
    # option processing
    input_filename = None
    for option, value in opts:
        if option in ('-i', '--input'):
            input_filename = value

    if not input_filename:
        print 'Paste sequence you want to cut:',
        seq = raw_input()
    else:
        seq = open(input_filename).read().strip()
    
    PrintPrimers(seq, START_HIS, INV_CAP_LINKER)
    
    
if __name__ == '__main__':
    main()
    
    