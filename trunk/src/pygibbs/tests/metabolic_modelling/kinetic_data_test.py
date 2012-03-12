#!/usr/bin/python

import unittest
import numpy as np

import itertools

from pygibbs.metabolic_modelling import kinetic_data


class UniformKineticDataTest(unittest.TestCase):
    
    def testBasic(self):
        kcat = 5
        km = 1e-5
        kdata = kinetic_data.UniformKineticData(kcat, km)
        
        for i in xrange(50):
            self.assertEqual(kdata.GetKcat(i), kcat)
        
        for i, j in itertools.product(range(10),range(10)):
            self.assertEqual(kdata.GetKm(i,j), km)
            
        kcats = np.matrix(np.ones((1, 10))) * kcat
        self.assertTrue((kcats == kdata.GetKcats(range(10))).all())

        kms = np.matrix(np.ones((10, 10))) * km
        self.assertTrue((kms == kdata.GetKm(range(10), range(10))).all())
        

def Suite():
    suites = (unittest.makeSuite(UniformKineticDataTest,'test'), )
    return unittest.TestSuite(suites)
    

if __name__ == '__main__':
    unittest.main()