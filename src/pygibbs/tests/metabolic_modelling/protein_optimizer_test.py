#!/usr/bin/python

import unittest
import numpy as np

from pygibbs.metabolic_modelling import kinetic_data
from pygibbs.metabolic_modelling import optimized_pathway
from pygibbs.metabolic_modelling import protein_optimizer
from pygibbs.metabolic_modelling import bounds

from pygibbs.tests.metabolic_modelling.fake_stoich_model import FakeStoichModel
from pygibbs.tests.metabolic_modelling.fake_thermo_data import FakeThermoData
from pygibbs.tests.metabolic_modelling.fake_thermo_data import FakeInfeasibleThermoData


class TestFixedVariableInjector(unittest.TestCase):
    
    def testOverlappingBounds(self):
        lb = np.matrix([[0,1,0,0,3]])
        ub = np.matrix([[5,1,5,5,3]])
        variable_i = np.where(lb != ub)
        
        injector = protein_optimizer.FixedVariableInjector(
            lb, ub, ub)
        
        variable_lb = lb[variable_i]
        variable_ub = ub[variable_i]
        self.assertTrue((variable_lb == injector.GetVariableLowerBounds()).all())
        self.assertTrue((variable_ub == injector.GetVariableUpperBounds()).all())
        self.assertTrue(
            (variable_ub == injector.GetVariableInitialConds()).all())
        
        expected_lb = injector(variable_lb)
        expected_ub = injector(variable_ub)
        
        self.assertTrue((expected_lb == lb).all())
        self.assertTrue((expected_ub == ub).all())
    
    def testNonoverlappingBounds(self):
        lb = np.matrix([[0,1,0,0,3]])
        ub = np.matrix([[5,2,5,5,4]])
        variable_i = np.where(lb != ub)
        
        injector = protein_optimizer.FixedVariableInjector(
            lb, ub, ub)
        
        variable_lb = lb[variable_i]
        variable_ub = ub[variable_i]
        self.assertTrue((variable_lb == injector.GetVariableLowerBounds()).all())
        self.assertTrue((variable_ub == injector.GetVariableUpperBounds()).all())
        self.assertTrue(
            (variable_ub == injector.GetVariableInitialConds()).all())
        
        expected_lb = injector(variable_lb)
        expected_ub = injector(variable_ub)
        
        self.assertTrue((expected_lb == lb).all())
        self.assertTrue((expected_ub == ub).all())
    
    def testInitialCondsDontMatchBounds(self):
        lb   = np.matrix([[0,1,0,0,3]])
        ub   = np.matrix([[5,1,5,5,3]])
        init = np.matrix([[-1,1,5,6,4]])
        variable_i = np.where(lb != ub)
        
        injector = protein_optimizer.FixedVariableInjector(
            lb, ub, init)
        
        variable_lb   = lb[variable_i]
        variable_ub   = ub[variable_i]
        expected_init = np.matrix([[0, 5, 5]])
        
        self.assertTrue((variable_lb == injector.GetVariableLowerBounds()).all())
        self.assertTrue((variable_ub == injector.GetVariableUpperBounds()).all())        
        self.assertTrue(
            (expected_init == injector.GetVariableInitialConds()).all())
        
        expected_lb = injector(variable_lb)
        expected_ub = injector(variable_ub)
        
        self.assertTrue((expected_lb == lb).all())
        self.assertTrue((expected_ub == ub).all())


class TestEnzymeLevelFunc(unittest.TestCase):
    
    def testApproximateDenom(self):
        RT = protein_optimizer.RT
        saturating_dg = -2*RT
        
        dGr_tags = np.matrix([[saturating_dg/2,
                               2*saturating_dg,
                               saturating_dg,
                               3*saturating_dg,
                               saturating_dg/3]])
        
        expected = np.matrix([[0.5, 1.0, 1.0, 1.0, 1.0/3.0]])
        actual   = protein_optimizer.EnzymeLevelFunc.ApproximateDenom(dGr_tags)
        diff     = np.abs(expected - actual) 
        
        self.assertTrue((diff < 1e-6).all())        
        

class DummyInjector(protein_optimizer.FixedVariableInjector):
    """Dummy injector does nothing."""
    
    def __init__(self):
        pass
    
    def __call__(self, x):
        return x
    

class TestMinusDG(unittest.TestCase):
    
    def testBasic(self):
        stoich_model = FakeStoichModel()
        thermo = FakeThermoData()

        dG0_r_prime = thermo.GetDGrTagZero_ForModel(stoich_model)
        S = stoich_model.GetStoichiometricMatrix()
        Ncompounds, _ = S.shape
        injector = DummyInjector()
        
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
        
        
class TestBoundDiffs(unittest.TestCase):
    
    def testBasic(self):
        lb = np.matrix([[0,1,0,0,3]])
        ub = np.matrix([[5,1,5,5,3]])
        bounder = protein_optimizer.BoundDiffs(lb, ub)

        self.assertTrue((bounder(ub) >= 0).all())
        self.assertTrue((bounder(lb) >= 0).all())
        
        within_bounds = (np.matrix([[1, 1, 1, 1, 3]]),
                         np.matrix([[2, 1, 2, 2, 3]]),
                         np.matrix([[1, 1, 2, 3, 3]]))
        for x in within_bounds:
            self.assertTrue((bounder(x) >= 0).all())
            
        out_of_bounds = (np.matrix([[-1, 1, 0, 5, 3]]),
                         np.matrix([[0, 2, 0, 5, 3]]),
                         np.matrix([[0, 1, 6, 5, 3]]),
                         np.matrix([[0, 1, 0, -3, 3]]),
                         np.matrix([[0, 1, 0, 2, 9]]),
                         np.matrix([[0, -7, 0, 2, 9]]))
        
        for x in out_of_bounds:
            self.assertFalse((bounder(x) >= 0).all())
        

class TestProteinOptimizer(unittest.TestCase):
    
    def MyBounds(self):
        my_lb, my_ub = 1e-8, 15e-3
        return bounds.Bounds(default_lb=my_lb,
                             default_ub=my_ub)
    
    def testDummyProblem(self):
        stoich_model = FakeStoichModel()
        thermo = FakeThermoData()
        kdata = kinetic_data.UniformKineticData()
        
        opt = protein_optimizer.ProteinOptimizer(
            stoich_model, thermo, kdata)
        res = opt.FindOptimum()
        self.assertEqual(optimized_pathway.OptimizationStatus.SUCCESSFUL,
                         res.status.status)
        
        result = res.opt_val
        self.assertAlmostEqual(0.048032, result, 3)
    
    def testDummyProblemDifferentBounds(self):
        stoich_model = FakeStoichModel()
        thermo = FakeThermoData()
        kdata = kinetic_data.UniformKineticData()
        
        b = self.MyBounds()    
        opt = protein_optimizer.ProteinOptimizer(
            stoich_model, thermo, kdata)
        res = opt.FindOptimum(concentration_bounds=b)
        self.assertEqual(optimized_pathway.OptimizationStatus.SUCCESSFUL,
                         res.status.status)
        
        result = res.opt_val
        self.assertAlmostEqual(0.041342, result, 3)
        
    def testInfeasibleDummyProblem(self):
        stoich_model = FakeStoichModel()
        thermo = FakeInfeasibleThermoData()
        kdata = kinetic_data.UniformKineticData()
        
        opt = protein_optimizer.ProteinOptimizer(
            stoich_model, thermo, kdata)
        res = opt.FindOptimum()
        self.assertEqual(optimized_pathway.OptimizationStatus.INFEASIBLE,
                         res.status.status)
        

def Suite():
    suites = (unittest.makeSuite(TestFixedVariableInjector,'test'),
              unittest.makeSuite(TestMinusDG,'test'),
              unittest.makeSuite(TestProteinOptimizer,'test'))
    return unittest.TestSuite(suites)
    

if __name__ == '__main__':
    unittest.main()