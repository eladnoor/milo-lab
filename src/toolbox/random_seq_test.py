#!/usr/bin/python

import unittest

from Bio.Alphabet import Alphabet
from Bio.Alphabet import IUPAC
from toolbox import random_seq


class MockAlphabet(Alphabet):
    
    def ContainsString(self, test_str):
        letter_set = set(self.letters).union(set(self.letters.upper()))
        for letter in test_str:
            if letter not in letter_set:
                return False
            
        return True


class OneLetterMockAlphabet(MockAlphabet):
    letters = 'a'


class TwoLetterMockAlphabet(MockAlphabet):
    letters = 'ab'
    
    
class MockDNAAlphabet(MockAlphabet):
    """Implemented purely for checking DNA matching using the base class."""
    letters = 'actg'


class TestRandomSeq(unittest.TestCase):
    
    def testWithOneLetterMockAlphabet(self):
        alphabet = OneLetterMockAlphabet()
        
        seq = random_seq.RandomSeq(0, alphabet=alphabet)
        self.assertEqual('', seq.tostring())
        
        seq = random_seq.RandomSeq(10, alphabet=alphabet)
        self.assertEqual('A'*10, seq.tostring())
        
        seq = random_seq.RandomSeq(100, alphabet=alphabet)
        self.assertEqual('A'*100, seq.tostring())

    def testWithTwoLetterMockAlphabet(self):
        alphabet = TwoLetterMockAlphabet()
        
        seq = random_seq.RandomSeq(0, alphabet=alphabet)
        self.assertEqual(0, len(seq.tostring()))
        self.assertEqual('', seq.tostring())
        
        seq = random_seq.RandomSeq(10, alphabet=alphabet)
        seq_str = seq.tostring()
        self.assertEqual(10, len(seq_str))
        self.assertTrue(alphabet.ContainsString(seq_str))
        
        seq = random_seq.RandomSeq(100, alphabet=alphabet)
        seq_str = seq.tostring()
        self.assertEqual(100, len(seq_str))
        self.assertTrue(alphabet.ContainsString(seq_str))
    
    def testWithDNAAlphabet(self):
        mock_alphabet = MockDNAAlphabet()
        alphabet = IUPAC.unambiguous_dna
        
        seq = random_seq.RandomSeq(0, alphabet=alphabet)
        self.assertEqual(0, len(seq.tostring()))
        self.assertEqual('', seq.tostring())
        
        seq = random_seq.RandomSeq(10, alphabet=alphabet)
        seq_str = seq.tostring()
        self.assertEqual(10, len(seq_str))
        self.assertTrue(mock_alphabet.ContainsString(seq_str))
        
        seq = random_seq.RandomSeq(100, alphabet=alphabet)
        seq_str = seq.tostring()
        self.assertEqual(100, len(seq_str))
        self.assertTrue(mock_alphabet.ContainsString(seq_str))
    
        
def Suite():
    return unittest.makeSuite(TestRandomSeq,'test')
    

if __name__ == '__main__':
    unittest.main()