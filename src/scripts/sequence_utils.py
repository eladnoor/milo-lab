
from Bio.Seq import MutableSeq


def MutableSeqFromFile(filename, alphabet):
    sequence_str = open(filename).read().strip()
    return MutableSeq(sequence_str.lower(), alphabet)
    

def DropStartCodon(seq):
    if seq.tostring().upper().startswith('ATG'):
        return seq[3:]
    return seq