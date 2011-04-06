#!/usr/bin/python

import random

from Bio.Alphabet import IUPAC
from Bio.Seq import Seq


class RandomSeq(Seq):
    
    @staticmethod
    def SampleLetters(letters, length):
        """Sample "length" letters from the "letters" set at random.
        
        Sampling with replacement since Python doesn't have a helper for it.
        
        Args:
            letters: the letters to choose.
            length: the number to choose.
            
        Returns:
            A string of the sampled letters in order.
        """
        chosen_letters = [None]*length
        for i in xrange(length):
            index = random.randint(0, len(letters)-1)
            chosen_letters[i] = letters[index]
        return ''.join(chosen_letters)
    
    def __init__(self, length, alphabet=IUPAC.unambiguous_dna):
        """Initialize a randomized sequence of the given length.
        
        Args:
            length: the sequence length.
            alphabet: the alphabet to choose from.
        """
        seq_str = self.SampleLetters(alphabet.letters, length)
        
        Seq.__init__(self, seq_str.upper(), alphabet)