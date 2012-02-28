#!/usr/bin/python

import unittest
import numpy as np

from pygibbs.metabolic_modelling import protein_optimizer
from pygibbs.metabolic_modelling import bounds

from pygibbs.tests.metabolic_modelling.fake_stoich_model import FakeStoichModel
from pygibbs.tests.metabolic_modelling.fake_thermo_data import FakeThermoData


class TestMinusDG(unittest.TestCase):
    
    def testBasic(self):
        stoich_model = FakeStoichModel()
        thermo = FakeThermoData()

        dG0_r_prime = thermo.GetDGrTagZero_ForModel(stoich_model)
        S = stoich_model.GetStoichiometricMatrix()
        Ncompounds, unused_Nreactions = S.shape
        injector = protein_optimizer.FixedVariableInjector(
            range(Ncompounds), [], [])
        
        minus_dg = protein_optimizer.MinusDG(S, dG0_r_prime, injector)
        
        # 1M concentrations should have no effect.
        x = np.ones(Ncompounds)
        transformed = minus_dg(x)
        expected_minus_dg = -np.matrix(dG0_r_prime)
        self.assertTrue((transformed == expected_minus_dg).all())
        
        # Biological concentrations should effect the dGrs.
        x = np.ones(Ncompounds) * 1e-3
        transformed = minus_dg(x)
        expected_minus_dg = -np.matrix(
            dG0_r_prime + protein_optimizer.RT * x * S)
        self.assertTrue((transformed == expected_minus_dg).all())
        

class TestProteinOptimizer(unittest.TestCase):
    
    def MyBounds(self):
        my_lb, my_ub = 1e-8, 15e-3
        return bounds.Bounds(default_lb=my_lb,
                             default_ub=my_ub)
    
    def testDummyProblem(self):
        stoich_model = FakeStoichModel()
        thermo = FakeThermoData()
        
        opt = protein_optimizer.ProteinOptimizer(stoich_model,
                                                 thermo)
        res = opt.FindOptimum()
        
        result = res.opt_val
        self.assertAlmostEqual(0.028924, result, 3)
    
    def testDummyProblemDifferentBounds(self):
        stoich_model = FakeStoichModel()
        thermo = FakeThermoData()
        
        b = self.MyBounds()    
        opt = protein_optimizer.ProteinOptimizer(stoich_model,
                                                 thermo)
        res = opt.FindOptimum(concentration_bounds=b)
        
        result = res.opt_val
        self.assertAlmostEqual(0.025950, result, 3)


def Suite():
    suites = (unittest.makeSuite(TestMinusDG,'test'),
              unittest.makeSuite(TestProteinOptimizer,'test'))
    return unittest.TestSuite(suites)
    

if __name__ == '__main__':
    unittest.main()