#!/usr/bin/python

import unittest
import numpy as np

from pygibbs.metabolic_modelling import mtdf_optimizer
from pygibbs.metabolic_modelling import bounds


class FakeStoichModel(object):
    
    def GetStoichiometricMatrix(self):
        return np.array([[-1,1,0],[0,-1,1]])
    
    def GetReactionIDs(self):
        return ['R1', 'R2']
    
    def GetCompoundIDs(self):
        return ['C1', 'C2', 'C3']
    
    def GetFluxes(self):
        return [1.0, 1.0]
    

class FakeThermoData(object):
    
    def GetDGrPrime(self, unused_s):
        return np.array([[11.1],
                         [-2.1]])
    

class TestMTDFOptimizer(unittest.TestCase):
    
    def MyBounds(self):
        my_lb, my_ub = 1e-8, 15e-3
        return bounds.Bounds(default_lb=my_lb,
                             default_ub=my_ub)
    
    def testDummyProblem(self):
        stoich_model = FakeStoichModel()
        thermo = FakeThermoData()

        S = stoich_model.GetStoichiometricMatrix()
        dg0r_primes = thermo.GetDGrPrime(S)
        
        opt = mtdf_optimizer.MTDFOptimizer(stoich_model, thermo)
        conc, mtdf = opt.FindMTDF()
        ln_conc = np.log(conc)
        transformed_dgr = dg0r_primes + mtdf_optimizer.RT * np.dot(S, ln_conc)
        expected_mtdf = -np.min(transformed_dgr)
        self.assertAlmostEqual(6.90989, mtdf, 3)
        self.assertAlmostEqual(expected_mtdf, mtdf, 3)
    
    def testDummyProblemDifferentBounds(self):
        stoich_model = FakeStoichModel()
        thermo = FakeThermoData()
        
        S = stoich_model.GetStoichiometricMatrix()
        dg0r_primes = thermo.GetDGrPrime(S)
        
        b = self.MyBounds()        
        opt = mtdf_optimizer.MTDFOptimizer(stoich_model, thermo)
        conc, mtdf = opt.FindMTDF(concentration_bounds=b)
        ln_conc = np.log(conc)
        transformed_dgr = dg0r_primes + mtdf_optimizer.RT * np.dot(S, ln_conc)
        expected_mtdf = -np.min(transformed_dgr)
        self.assertAlmostEqual(13.117132, mtdf, 3)
        self.assertAlmostEqual(expected_mtdf, mtdf, 3)


def Suite():
    return unittest.makeSuite(TestMTDFOptimizer,'test')
    

if __name__ == '__main__':
    unittest.main()