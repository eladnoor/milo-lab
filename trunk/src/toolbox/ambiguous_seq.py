#!/usr/bin/python

import itertools

from Bio.Alphabet import DNAAlphabet
from Bio.Alphabet.IUPAC import IUPACAmbiguousDNA
from Bio.Seq import Seq


class AmbigousDNASeq(Seq):
    """Class containing ambiguous DNA sequences."""
    
    AMBIGUOUS_BASE_MAP = {'R': ('G', 'A'),  # Purines
                          'Y': ('T', 'C'),  # Pyrimidine
                          'K': ('G', 'T'),  # Keto
                          'M': ('A', 'C'),  # Amino
                          'S': ('G', 'C'),
                          'W': ('A', 'T'),
                          'B': ('G', 'T', 'C'),
                          'D': ('G', 'A', 'T'),
                          'H': ('A', 'C', 'T'),
                          'V': ('G', 'C', 'A'),
                          'N': ('A', 'G', 'C', 'T'),
                          }
    
    def __init__(self, seq_str):
        Seq.__init__(self, seq_str.upper(), IUPACAmbiguousDNA())
        
        self._seq_list = list(self.tostring())
    
    def AmbiguousPositions(self):
        """Return which positions have ambiguous bases."""
        return [i for i, base in enumerate(self._seq_list)
                if base in self.AMBIGUOUS_BASE_MAP]
        
    def IsConcrete(self):
        """Returns True if this sequence is not ambiguous."""
        return not self.AmbiguousPositions()
    
    def GetConcreteSeq(self):
        """Returns a Seq object if this is a concrete sequence."""
        if not self.IsConcrete():
            raise TypeError('This sequence %s is not concrete' % self.tostring())
        
        return Seq(self.tostring(), DNAAlphabet())
    
    def AllConcreteSeqs(self):
        """An iterator over all the concrete sequences denoted by this
           ambiguous one.
        
        Yields:
            Each possible concrete sequence as a Bio.Seq.Seq object.
        """
        ambiguous_positions = self.AmbiguousPositions()
        
        ambiguous_position_ranges = []
        for pos in ambiguous_positions:
            ambig_base = self._seq_list[pos]
            substitutions = self.AMBIGUOUS_BASE_MAP[ambig_base]
            r = range(len(substitutions))
            ambiguous_position_ranges.append(range(len(substitutions)))
        
        for concrete_indices in itertools.product(*ambiguous_position_ranges):
            new_seq_list = list(self._seq_list)
            
            for ambig_index, concrete_index in enumerate(concrete_indices):
                pos = ambiguous_positions[ambig_index]
                ambig_base = self._seq_list[pos]
                substitution = self.AMBIGUOUS_BASE_MAP[ambig_base][concrete_index]
                
                new_seq_list[pos] = substitution
            
            yield Seq(''.join(new_seq_list), DNAAlphabet())
                
                
