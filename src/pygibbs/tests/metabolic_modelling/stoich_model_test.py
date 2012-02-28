#!/usr/bin/python

import unittest
import numpy as np

from pygibbs.metabolic_modelling import stoich_model


class TestStoichiometricModel(unittest.TestCase):
    
    def testWrongNumberOfReactions(self):
        # Fake matrix with 5 reactions and 5 compounds 
        S = np.ones((5,5))
        compound_ids = range(8)
        reaction_ids = range(5)
        
        self.assertRaises(ValueError, stoich_model.StoichiometricModel,
                          S, reaction_ids, compound_ids)

    def testWrongNumberOfCompounds(self):
        # Fake matrix with 5 reactions and 5 compounds 
        S = np.ones((5,5))
        compound_ids = range(5)
        reaction_ids = range(11)
        
        self.assertRaises(ValueError, stoich_model.StoichiometricModel,
                          S, reaction_ids, compound_ids)        

    def testBasic(self):
        # Fake matrix with 5 reactions and 5 compounds 
        S = np.ones((5,5))
        compound_ids = range(5)
        reaction_ids = range(5)
        
        model = stoich_model.StoichiometricModel(
            S, reaction_ids, compound_ids)
        equality_m = (S == model.GetStoichiometricMatrix())
        self.assertTrue(equality_m.all())
        self.assertEqual(compound_ids, model.GetCompoundIDs())
        self.assertEqual(reaction_ids, model.GetReactionIDs())
        
        fluxes = model.GetFluxes()
        self.assertEqual(len(reaction_ids), fluxes.size)
        expected_fluxes = np.ones((1, len(reaction_ids)))
        self.assertTrue((expected_fluxes == fluxes).all())
        

def Suite():
    suites = (unittest.makeSuite(TestStoichiometricModel,'test'),)
    return unittest.TestSuite(suites)
    

if __name__ == '__main__':
    unittest.main()