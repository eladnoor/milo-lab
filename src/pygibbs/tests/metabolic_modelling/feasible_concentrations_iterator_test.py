#!/usr/bin/python

import unittest
import numpy as np

from pygibbs.metabolic_modelling import feasible_concentrations_iterator
from pygibbs.metabolic_modelling import bounds

from pygibbs.tests.metabolic_modelling.fake_stoich_model import FakeStoichModel
from pygibbs.tests.metabolic_modelling.fake_thermo_data import FakeThermoData
    

class TestFeasibleConcentrationsIterator(unittest.TestCase):
    
    def NormalBounds(self):
        my_lb, my_ub = 1e-6, 1e-2
        return bounds.Bounds(default_lb=my_lb,
                             default_ub=my_ub)
    
    def WideBounds(self):
        my_lb, my_ub = 1e-8, 20e-3
        return bounds.Bounds(default_lb=my_lb,
                             default_ub=my_ub)
    
    def testDummyProblem(self):
        stoich_model = FakeStoichModel()
        thermo = FakeThermoData()
        b = self.NormalBounds()
        
        iter = feasible_concentrations_iterator.FeasibleConcentrationsIterator(
            stoich_model, thermo, b)
        for concs in iter:
            self.assertTrue(iter.Feasible(concs))
            
    def testDummyProblemWideBounds(self):
        stoich_model = FakeStoichModel()
        thermo = FakeThermoData()
        b = self.WideBounds()
        
        iter = feasible_concentrations_iterator.FeasibleConcentrationsIterator(
            stoich_model, thermo, b)
        for concs in iter:
            self.assertTrue(iter.Feasible(concs))


def Suite():
    return unittest.makeSuite(TestFeasibleConcentrationsIterator,'test')
    

if __name__ == '__main__':
    unittest.main()