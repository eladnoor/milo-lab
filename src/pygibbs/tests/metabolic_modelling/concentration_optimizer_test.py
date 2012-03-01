#!/usr/bin/python

import unittest
import numpy as np

from pygibbs.metabolic_modelling import concentration_optimizer
from pygibbs.metabolic_modelling import bounds

from pygibbs.tests.metabolic_modelling.fake_stoich_model import FakeStoichModel
from pygibbs.tests.metabolic_modelling.fake_thermo_data import FakeThermoData
    

class TestConcentrationOptimizer(unittest.TestCase):
    
    def MyBounds(self):
        my_lb, my_ub = 1e-8, 15e-3
        return bounds.Bounds(default_lb=my_lb,
                             default_ub=my_ub)
    
    def testDummyProblem(self):
        stoich_model = FakeStoichModel()
        thermo = FakeThermoData()
        
        opt = concentration_optimizer.ConcentrationOptimizer(stoich_model, thermo)
        
        res = opt.MinimizeConcentration(0)
        val = res.opt_val
        self.assertAlmostEqual(8.82430e-5, val, 3)
        
        res = opt.MinimizeConcentration(2)
        val = res.opt_val
        
        # Makes sense - last concentration can be as low as possible.
        self.assertAlmostEqual(1.0e-6, val, 3)
    
    def testDummyProblemDifferentBounds(self):
        stoich_model = FakeStoichModel()
        thermo = FakeThermoData()
        
        b = self.MyBounds()        
        opt = concentration_optimizer.ConcentrationOptimizer(stoich_model, thermo)
        res = opt.MinimizeConcentration(0, concentration_bounds=b)
    
        val = res.opt_val
        self.assertAlmostEqual(8.82478630613e-7, val, 3)
    
        res = opt.MinimizeConcentration(2)
        val = res.opt_val
        
        # Makes sense - last concentration can be as low as possible.
        self.assertAlmostEqual(1.0e-8, val, 3)

def Suite():
    return unittest.makeSuite(TestConcentrationOptimizer,'test')
    

if __name__ == '__main__':
    unittest.main()