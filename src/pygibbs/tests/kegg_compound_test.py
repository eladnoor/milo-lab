#!/usr/bin/python

import unittest

from pygibbs import kegg_compound
from pygibbs import pseudoisomer

class TestKeggCompound(unittest.TestCase):
    
    def setUp(self):
        self.cid = 1323
        self.str_cid = str(self.cid)
        self.name = 'Glyxolaterium-monohysterianase'
        self.names = ['Glyxolaterium-monohysterianase',
                      'GMH']
        self.mass = '60.99'
        self.formula = 'H2CO3'
        self.inchi = 'InChI=1S/C2H2O2/c3-1-2-4/h1-2H' # inchi for glyoxal
        
        self.test_compound = kegg_compound.Compound(
            cid=self.cid, name=self.name, all_names=self.names,
            mass=self.mass, formula=self.formula, inchi=self.inchi)
        
    def testProperties(self):
        self.assertEqual(self.cid, self.test_compound.cid)
        self.assertEqual(self.name, self.test_compound.name)
        self.assertEqual(self.names, self.test_compound.all_names)
        self.assertEqual(self.mass, self.test_compound.mass)
        self.assertEqual(self.formula, self.test_compound.formula)
        self.assertEqual(self.inchi, self.test_compound.inchi)
        self.assertTrue(self.test_compound.from_kegg)
        
        kegg_link = self.test_compound.kegg_link
        self.assertNotEqual(-1, kegg_link.find(self.str_cid))
        
    def testSimpleGetters(self):
        molecule = self.test_compound.GetMolecule()
        self.assertTrue(molecule is not None)
        
        smiles = self.test_compound.get_smiles()
        self.assertTrue(smiles is not None)
        
        self.assertEqual(self.inchi, self.test_compound.get_inchi())
        self.assertNotEqual(-1, self.test_compound.get_string_cid().find(self.str_cid))
        
        # We call these methods to make sure they don't error.
        self.test_compound.get_atom_bag()
        self.test_compound.get_nH_and_charge()
        self.test_compound.get_num_electrons()
        
    def testAddThermodynamicData(self):
        pmap = pseudoisomer.PseudoisomerMap(nH=8, z=0, nMg=0, dG0=18.8)
        source_string = 'Noor, Bar-Even 2011'
        self.test_compound.AddThermodynamicData(pmap, source_string)
        self.assertEqual(pmap, self.test_compound.pmap)
        self.assertEqual(source_string, self.test_compound.pmap_source)
    
    def testSetThermodynamicError(self):
        error = 'Group contribution failed to solve your problems.'
        self.test_compound.SetThermodynamicError(error)
        self.assertEqual(error, self.test_compound.pmap_error)
    
    def testJSONDict(self):
        json_dict = self.test_compound.ToJSONDict()
        
        key_set = set(json_dict.keys())
        expected_keys = set(['CID', 'mass', 'formula',
                             'names', 'InChI', 'num_electrons'])
        self.assertEqual(expected_keys, key_set)
        
        self.assertEqual(self.test_compound.get_string_cid(), json_dict['CID'])
        self.assertEqual(self.mass, json_dict['mass'])
        self.assertEqual(self.formula, json_dict['formula'])
        self.assertEqual(self.names, json_dict['names'])
        self.assertEqual(self.inchi, json_dict['InChI'])
        self.assertTrue(json_dict['num_electrons'] is not None)
    
    def testToDBRow(self):
        expected_row = [self.cid, self.name, ';'.join(self.names),
                        self.mass, self.formula, self.inchi,
                        self.test_compound.get_num_electrons(),
                        self.test_compound.from_kegg,
                        self.test_compound.pubchem_id,
                        self.test_compound.cas]
        db_row = self.test_compound.ToDBRow()
        self.assertEqual(expected_row, db_row)
        

def Suite():
    return unittest.makeSuite(TestKeggCompound,'test')
    

if __name__ == '__main__':
    unittest.main()