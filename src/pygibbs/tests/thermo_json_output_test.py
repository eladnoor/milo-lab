#!/usr/bin/python

import csv
import unittest

from StringIO import StringIO
from pygibbs.thermodynamics import PsuedoisomerTableThermodynamics
from toolbox.database import SqliteDatabase


CSV_DATA = """"smiles","cid","compound name","dG0","dH0","z","nH","nMg","use for","ref","remark"
,1,"H2O",-237.19,-285.83,0,2,0,"test","Alberty (2003)",
,2,"ATP",-2768.1,-3619.21,-4,12,0,"train","Alberty (2003)",
,2,"ATP",-3258.7,-4063.31,-2,12,1,"test","Alberty (2006)","with Mg"
,7,"O2(aq)",16.4,-11.7,0,0,0,"test","Alberty (2003)",
,7,"O2(g)",0,0,0,0,0,"skip","Alberty (2003)",
"""

PUBLIC_DB_FNAME = '../data/public_data.sqlite'

class TestThermoJsonOutput(unittest.TestCase):
    
    def setUp(self):
        fake_csv_file = StringIO(CSV_DATA)
        csv_reader = csv.DictReader(fake_csv_file)
        self.fake_thermo_csv = PsuedoisomerTableThermodynamics()
        self.fake_thermo_csv = PsuedoisomerTableThermodynamics._FromDictReader(
                                    csv_reader, self.fake_thermo_csv,
                                    warn_for_conflicting_refs=False)
        
        db = SqliteDatabase(PUBLIC_DB_FNAME)
        db_reader = db.DictReader('fake_pseudoisomers')
        self.fake_thermo_db = PsuedoisomerTableThermodynamics()
        self.fake_thermo_db = PsuedoisomerTableThermodynamics._FromDictReader(
                                    db_reader, self.fake_thermo_db,
                                    warn_for_conflicting_refs=False)

    def testGetJsonDictionary(self):
        json_list = [self.fake_thermo_csv.GetJSONDictionary(),
                     self.fake_thermo_db.GetJSONDictionary()]
        
        for json_data in json_list:
            self.assertEqual(json_data[0]['cid'], 1)
            self.assertEqual(json_data[0]['source'], 'Alberty (2003)')
            self.assertEqual(json_data[0]['inchi'], 'InChI=1S/H2O/h1H2')
            self.assertEqual(json_data[0]['num_electrons'], 10)
            self.assertAlmostEqual(json_data[0]['species'][0]['dG0_f'], -237.19, places=5)
    
            self.assertEqual(json_data[1]['cid'], 2)
            self.assertEqual(json_data[1]['source'], 'Alberty (2003)')
            self.assertEqual(json_data[1]['inchi'], 'InChI=1S/C10H16N5O13P3/c11-8-5-9(13-2-12-8)15(3-14-5)10-7(17)6(16)4(26-10)1-25-30(21,22)28-31(23,24)27-29(18,19)20/h2-4,6-7,10,16-17H,1H2,(H,21,22)(H,23,24)(H2,11,12,13)(H2,18,19,20)/t4-,6-,7-,10-/m1/s1')
            self.assertEqual(json_data[1]['num_electrons'], 260)
            self.assertEqual(json_data[1]['species'][0]['nH'], 12)
            self.assertEqual(json_data[1]['species'][0]['nMg'], 0)
            self.assertEqual(json_data[1]['species'][0]['z'], -4)
            self.assertAlmostEqual(json_data[1]['species'][0]['dG0_f'], -2768.1, places=5)
            self.assertEqual(json_data[1]['species'][1]['nH'], 12)
            self.assertEqual(json_data[1]['species'][1]['nMg'], 1)
            self.assertEqual(json_data[1]['species'][1]['z'], -2)
            self.assertAlmostEqual(json_data[1]['species'][1]['dG0_f'], -3258.7, places=5)
    
            self.assertEqual(json_data[2]['cid'], 7)
            self.assertEqual(json_data[2]['source'], 'Alberty (2003)')
            self.assertEqual(json_data[2]['inchi'], 'InChI=1S/O2/c1-2')
            self.assertEqual(json_data[2]['num_electrons'], 16)
            self.assertAlmostEqual(json_data[2]['species'][0]['dG0_f'], 16.4, places=5)

def Suite():
    return unittest.makeSuite(TestThermoJsonOutput, 'test')
    

if __name__ == '__main__':
    unittest.main()