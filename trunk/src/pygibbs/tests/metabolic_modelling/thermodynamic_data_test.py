#!/usr/bin/python

import unittest
import numpy as np

from pygibbs.metabolic_modelling import thermodynamic_data


class FakeStoichModel(object):
    
    def GetStoichiometricMatrix(self):
        return np.array([[-1,1,0],[0,-1,1]])
    
    def GetReactionIDs(self):
        return ['R1', 'R2']
    
    def GetCompoundIDs(self):
        return ['C1', 'C2', 'C3']
    
    def GetFluxes(self):
        return [1.0, 1.0]
    

class TestFormationBasedThermoData(unittest.TestCase):
    
    def testBasic(self):
        formation_energies = {'C1':1,
                              'C2':-10,
                              'C3':395.120}
        
        thermo_data = thermodynamic_data.FormationBasedThermoData(formation_energies)
        for id, formation in formation_energies.iteritems():
            self.assertEquals(formation, thermo_data.GetDGfTagZero_ForID(id))
        
        self.assertTrue(np.isnan(thermo_data.GetDGfTagZero_ForID('FAKE_ID')))
        
        # Fetch compounds that have data.
        compound_ids = formation_energies.keys()
        expected_formation = [formation_energies[k] for k in compound_ids]
        actual_formation = thermo_data.GetDGfTagZero_ForIDs(compound_ids)
        self.assertEqual(expected_formation, list(actual_formation))
        
        # Fetch mixture of compounds, one that has no data.
        compound_ids.append('FAKE_ID')
        actual_formation = list(thermo_data.GetDGfTagZero_ForIDs(compound_ids))
        
        self.assertEqual(expected_formation, actual_formation[:-1])
        self.assertTrue(np.isnan(actual_formation[-1]))

        # Check that you can't get reaction energies by ID.
        self.assertRaises(NotImplementedError,
                          thermo_data.GetDGrTagZero_ForID, 
                          'FAKE_ID')
        
        # Check that we can get thermo data from the model.
        model = FakeStoichModel()
        S = model.GetStoichiometricMatrix()
        formation_energies = thermo_data.GetDGfTagZero_ForModel(model)
        expected_reaction_energies = list(np.dot(S, formation_energies))
        reaction_energies = list(thermo_data.GetDGrTagZero_ForModel(model))
        self.assertEqual(expected_reaction_energies, reaction_energies)


class TestReactionBasedThermoData(unittest.TestCase):
    
    def testBasic(self):
        reaction_energies = {'R1':1,
                             'R2':-10,
                             'R3':395.120}
        
        thermo_data = thermodynamic_data.ReactionBasedThermoData(reaction_energies)
        for id, r_energy in reaction_energies.iteritems():
            self.assertEquals(r_energy, thermo_data.GetDGrTagZero_ForID(id))
        
        self.assertTrue(np.isnan(thermo_data.GetDGrTagZero_ForID('FAKE_ID')))        
        
        # Fetch compounds that have data.
        reaction_ids = reaction_energies.keys()
        expected_energies = [reaction_energies[k] for k in reaction_ids]
        actual_energies = thermo_data.GetDGrTagZero_ForIDs(reaction_ids)
        self.assertEqual(expected_energies, list(actual_energies))
                
        # Fetch mixture of compounds, one that has no data.
        reaction_ids.append('FAKE_ID')
        actual_energies = list(thermo_data.GetDGrTagZero_ForIDs(reaction_ids))
        
        self.assertEqual(expected_energies, actual_energies[:-1])
        self.assertTrue(np.isnan(actual_energies[-1]))

        # Check that you can't get compound energies by ID.
        self.assertRaises(NotImplementedError,
                          thermo_data.GetDGfTagZero_ForID, 
                          'FAKE_ID')
        
        # Check that we can get thermo data from the model.
        model = FakeStoichModel()
        reaction_ids = model.GetReactionIDs()
        expected_energies = thermo_data.GetDGrTagZero_ForIDs(reaction_ids)
        actual_energies = thermo_data.GetDGrTagZero_ForModel(model)
        self.assertEqual(list(expected_energies), list(actual_energies))
        

def Suite():
    suites = (unittest.makeSuite(TestFormationBasedThermoData, 'test'),
              unittest.makeSuite(TestReactionBasedThermoData, 'test'))
    return unittest.TestSuite(suites)
    

if __name__ == '__main__':
    unittest.main()