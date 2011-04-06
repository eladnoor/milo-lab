#!/usr/bin/python

import unittest

from toolbox import ambiguous_seq


class TestAmbigousDNASeq(unittest.TestCase):
    
    def assertSeqsEqual(self, expected_seqs, actual_seqs):
        expected_set = set([str(s) for s in expected_seqs])
        actual_set = set([str(s) for s in actual_seqs])
        self.assertEqual(expected_set, actual_set)
    
    def testAmbiguousPositions(self):
        
        seq = ambiguous_seq.AmbigousDNASeq('ACTGGGA')
        self.assertEqual([], seq.AmbiguousPositions())
        self.assertTrue(seq.IsConcrete())
        self.assertEqual('ACTGGGA', seq.GetConcreteSeq().tostring())
        
        seq = ambiguous_seq.AmbigousDNASeq('ACNGGVA')
        self.assertEqual([2,5], seq.AmbiguousPositions())
        
        seq = ambiguous_seq.AmbigousDNASeq('ACNYGVACGTNNND')
        self.assertEqual([2, 3, 5, 10, 11, 12, 13], seq.AmbiguousPositions())
    
    def testAllConcreteSeqsSimple(self):
        ambig_seq = ambiguous_seq.AmbigousDNASeq('ACTGGGA')
        concrete_seq = ambig_seq.GetConcreteSeq()
        self.assertSeqsEqual([concrete_seq], ambig_seq.AllConcreteSeqs())
        
    def testAllConcreteSeqsWithAmbigBases(self):
        ambig_seq = ambiguous_seq.AmbigousDNASeq('ACNGGVA')
        expected_seqs = ['ACAGGGA', 'ACCGGGA', 'ACTGGGA', 'ACGGGGA',
                         'ACAGGCA', 'ACCGGCA', 'ACTGGCA', 'ACGGGCA',
                         'ACAGGAA', 'ACCGGAA', 'ACTGGAA', 'ACGGGAA']
        self.assertSeqsEqual(expected_seqs, ambig_seq.AllConcreteSeqs())


def Suite():
    return unittest.makeSuite(TestAmbigousDNASeq,'test')
    

if __name__ == '__main__':
    unittest.main()