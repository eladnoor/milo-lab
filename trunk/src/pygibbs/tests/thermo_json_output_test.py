#!/usr/bin/python

import csv
import unittest

from StringIO import StringIO
from pygibbs.thermodynamics import PsuedoisomerTableThermodynamics
from toolbox.database import SqliteDatabase


CSV_DATA = """"smiles","cid","compound name","dG0","dH0","z","nH","nMg","use for","ref","remark"
,1,"H2O",-237.19,-285.83,0,2,0,"test","Alberty (2003)",
"C([C@@H]1[C@H]([C@H]([C@H](n2cnc3c(N)ncnc23)O1)O)O)OP(=O)([O-])OP(=O)([O-])OP(=O)([O-])[O-]",2,"ATP",-2768.1,-3619.21,-4,12,0,"train","Alberty (2003)",
"C([C@@H]1[C@H]([C@H]([C@H](n2cnc3c(N)ncnc23)O1)O)O)OP(=O)([O-])OP(=O)([O-])OP(=O)([O-])O",2,"ATP",-2811.48,-3612.91,-3,13,0,"test","Alberty (2003)",
"C([C@@H]1[C@H]([C@H]([C@H](n2cnc3c(N)[nH+]cnc23)O1)O)O)OP(=O)([O-])OP(=O)([O-])OP(=O)([O-])O",2,"ATP",-2838.18,-3627.91,-2,14,0,"test","Alberty (2003)",
"C([C@@H]1[C@H]([C@H]([C@H](n2cnc3c(N)ncnc23)O1)O)O)OP(=O)([O-])OP(=O)([O-])OP(=O)([O-])[O-].[Mg+2]",2,"ATP",-3258.7,-4063.31,-2,12,1,"test","Alberty (2006)","with Mg"
"C([C@@H]1[C@H]([C@H]([C@H](n2cnc3c(N)[nH+]cnc23)O1)O)O)OP(=O)([O-])OP(=O)([O-])OP(=O)([O-])[O-].[Mg+2]",2,"ATP",-3287.5,-4063.01,-1,13,1,"test","Alberty (2006)","with Mg"
"C([C@@H]1[C@H]([C@H]([C@H](n2cnc3c(N)ncnc23)O1)O)O)OP(=O)([O-])OP(=O)([O-])OP(=O)([O-])[O-].[Mg+2].[Mg+2]",2,"ATP",-3729.3,-4519.51,0,12,2,"test","Alberty (2006)","with Mg"
,7,"O2(aq)",16.4,-11.7,0,0,0,"test","Alberty (2003)",
,7,"O2(g)",0,0,0,0,0,"skip","Alberty (2003)",
"""


class TestThermoJsonOutput(unittest.TestCase):
    
    def setUp(self):
        fake_csv_file = StringIO(CSV_DATA)
        reader = csv.DictReader(fake_csv_file)
        self.fake_thermo = PsuedoisomerTableThermodynamics._FromDictReader(reader)

        # TODO(elad): test using a database too.

    def testGetJsonDictionary(self):
        json_data = self.fake_thermo.get_json_dictionary()
        
        # TODO(elad): actually test data. 
        
        self.assertTrue(self.fake_thermo.cid2PseudoisomerMap(1) is not None)


def Suite():
    return unittest.makeSuite(TestThermoJsonOutput, 'test')
    

if __name__ == '__main__':
    unittest.main()