#!/usr/bin/python

import unittest
import numpy as np

from pygibbs.metabolic_modelling import mtdf_optimizer
from pygibbs.metabolic_modelling import bounds

from pygibbs.tests.metabolic_modelling.fake_stoich_model import FakeStoichModel
from pygibbs.tests.metabolic_modelling.fake_thermo_data import FakeThermoData
    

class TestMTDFOptimizer(unittest.TestCase):
    
    def MyBounds(self):
        my_lb, my_ub = 1e-8, 15e-3
        return bounds.Bounds(default_lb=my_lb,
                             default_ub=my_ub)
    
    def testDummyProblem(self):
        stoich_model = FakeStoichModel()
        thermo = FakeThermoData()
        
        opt = mtdf_optimizer.MTDFOptimizer(stoich_model, thermo)
        res = opt.FindMTDF()
        
        transformed_dgr = res.dGr_tag
        expected_mtdf = -np.min(transformed_dgr)
        mtdf = res.opt_val
        
        self.assertAlmostEqual(6.90989, mtdf, 3)
        self.assertAlmostEqual(expected_mtdf, mtdf, 3)
    
    def testDummyProblemDifferentBounds(self):
        stoich_model = FakeStoichModel()
        thermo = FakeThermoData()
        
        b = self.MyBounds()
        opt = mtdf_optimizer.MTDFOptimizer(stoich_model, thermo)
        res = opt.FindMTDF(concentration_bounds=b)
        
        transformed_dgr = res.dGr_tag
        expected_mtdf = -np.min(transformed_dgr)
        mtdf = res.opt_val
        
        self.assertAlmostEqual(13.117132, mtdf, 3)
        self.assertAlmostEqual(expected_mtdf, mtdf, 3)


def Suite():
    return unittest.makeSuite(TestMTDFOptimizer,'test')
    

if __name__ == '__main__':
    unittest.main()