#!/usr/bin/python

import unittest

from pygibbs import kegg_enzyme


class TestKeggEnzyme(unittest.TestCase):
    
    def setUp(self):
        self.ec = '1.1.1.1'
        self.title = 'Transketolase QMFD reaction KDPG:LMFQ!'
        self.names = ['Ketolase D',
                      'Transketolase D']
        self.reactions = [123,234,1212]
        
        self.test_enzyme = kegg_enzyme.Enzyme(
            self.ec, title=self.title, names=self.names,
            reactions=self.reactions)
    
    def testProcessEC(self):
        expected_ec = '1.1.1.2'
        test_ec_strings = ('EC 1.1.1.2',
                           '1.1.1.2',
                           '  1.1.1.2',
                           '1.1.1.2  ',
                           '  1.1.1.2  ')
                
        for test_ec_str in test_ec_strings:
            self.assertEqual(expected_ec,
                             self.test_enzyme.ProcessEC(test_ec_str))
    
    def testGetCompoundIds(self):
        test_pairs = (('alcohol [CPD:C00069]; NAD+ [CPD:C00003]',
                       ['C00069', 'C00003']),
                      ('L-homoserine [CPD:C00263]; NAD+ [CPD:C00003]; NADP+ [CPD:C00006]',
                       ['C00263', 'C00003', 'C00006']),
                      ('UDP-glucuronate [CPD:C00167];  NADH [CPD:C00004];  H+ [CPD:C00080]',
                       ['C00167', 'C00004', 'C00080']),
                      ('UDP-glucuronate;  NADH ;  H+',
                       []))
        
        for test_str, expected_cids in test_pairs:
            res = self.test_enzyme.GetCompoundIds(test_str)
            self.assertEqual(expected_cids, res)
    
    def testGetStringRID(self):
        test_pairs = ((1234, 'R01234'),
                      (1, 'R00001'),
                      (51, 'R00051'))
        
        for int_id, str_id in test_pairs:
            self.assertEqual(str_id, self.test_enzyme.GetStringRID(int_id))
        
    def testProperties(self):
        self.assertEqual(self.ec, self.test_enzyme.ec)
        self.assertEqual(self.title, self.test_enzyme.title)
        self.assertEqual(self.names, self.test_enzyme.names)
        self.assertEqual(self.reactions, self.test_enzyme.reactions)
        
        kegg_link = self.test_enzyme.kegg_link
        self.assertNotEqual(-1, kegg_link.find(self.ec))
    
    def testJSONDict(self):
        json_dict = self.test_enzyme.ToJSONDict()
        
        key_set = set(json_dict.keys())
        expected_keys = set(['EC', 'title', 'names', 'reaction_ids',
                             'substrates', 'products', 'cofactors',
                             'organisms', 'orthology', 'genes'])
        self.assertEqual(expected_keys, key_set)
        
        self.assertEqual(self.ec, json_dict['EC'])
        self.assertEqual(self.title, json_dict['title'])
        self.assertEqual(self.names, json_dict['names'])
        
        expected_rids = map(self.test_enzyme.GetStringRID, self.reactions)
        self.assertEqual(expected_rids, json_dict['reaction_ids'])


def Suite():
    return unittest.makeSuite(TestKeggEnzyme, 'test')
    

if __name__ == '__main__':
    unittest.main()