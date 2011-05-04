#!/usr/bin/python

import StringIO
import unittest

from pygibbs import pathway
from pygibbs.kegg_parser import ParsedKeggFile


class TestPathwayConditions(unittest.TestCase):
    
    def testFromString(self):
        test_strings = ('pH=7.0,I=0.1,T=300',
                        'media=None,pH=7.0,I=0.1,T=300\t'
                        'media=glucose,pH=7.0,I=0.1,T=300',
                        'media=glucose,pH=7.0,I=0.1,T=300,c0=0.0001',
                        'media=None,pH=7.0,I=0.1,T=300,c0=0.0001\t'
                        'media=glucose,pH=7.0,I=0.1,T=300,c0=0.0001')
        
        for t_str in test_strings:
            conds = pathway.PathwayConditions.FromString(t_str)
            self.assertTrue(len(conds) > 0)


TEST_FILE_CONTENTS = """
ENTRY       L-CALVIN_NAD
SKIP        TRUE
NAME        rPP cycle (NAD)
TYPE        PROFILE
CONDITIONS  media=None,pH=7.0,I=0.1,T=300,c0=0.0001
            media=glucose,pH=7.0,I=0.1,T=300,c0=0.0001
REACTION    R00024  C01182 + C00011 + C00001 -> 2 C00197
            R01512  C00002 + C00197 -> C00008 + C00236
            R01061  C00236 + C00004 + C00080 -> C00118 + C00009 + C00003
            R01015  C00118 -> C00111
            R01068  C00111 + C00118 -> C00354
            R00762  C00354 + C00001 -> C00085 + C00009
            R01067  C00085 + C00118 -> C00279 + C00231
            R01829  C00111 + C00279 -> C00447
            R01845  C00447 + C00001 -> C05382 + C00009
            R01641  C05382 + C00118 -> C00117 + C00231
            R01056  C00117 -> C00199
            R01529  C00231 -> C00199 
            R01523  C00002 + C00199 -> C00008 + C01182
///
ENTRY       M00093
SKIP        TRUE
NAME        Uridine monophosphate biosynthesis
TYPE        MARGIN
CONDITIONS  pH=7.0,I=0.1,T=300
C_MID       1e-3
MODULE      M00093
MAP_CID     C00201 C00002
            C00454 C00008
///
ENTRY       C-GAPDH
SKIP        TRUE
NAME        GAP to 3PG
TYPE        CONTOUR
REACTION    C00118 + C00009 + C00008 + C00003 -> C00002 + C00004 + C00197 + C00080
PH          0.0 1.0 2.0 3.0 4.0 5.0 5.5 6.0 6.5 7.0 7.5 8.0 8.5 9.0 9.5 10.0
I           0.0 0.05 0.1 0.15 0.2 0.25 0.3 0.35 0.4
///
ENTRY       P-RuBP
SKIP        TRUE
NAME        pKa of Ribulose-1,5P
TYPE        PROTONATION
PH          5.0 5.2 5.4 5.6 5.8 6.0 6.2 6.4 6.6 6.8 7.0 7.2 7.4 7.6 7.9 8.0 8.2 8.4 8.6 8.8 9.0
COMPOUND    C01182
///
"""

class TestPathwayData(unittest.TestCase):
    
    def setUp(self):
        self.fake_file = StringIO.StringIO(TEST_FILE_CONTENTS)
    
    def testFromFieldMap(self):
        parsed_kegg = ParsedKeggFile._FromKeggFileHandle(self.fake_file)
        
        for field_map in parsed_kegg.values():
            data = pathway.PathwayData.FromFieldMap(field_map)
            self.assertNotEqual(None, data)
            

def Suite():
    suites = (unittest.makeSuite(TestPathwayConditions,'test'),
              unittest.makeSuite(TestPathwayData,'test'))
    return unittest.TestSuite(suites)

if __name__ == '__main__':
    unittest.main()