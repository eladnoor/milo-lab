#!/usr/bin/python

import unittest
import numpy as np

from pygibbs.metabolic_modelling import bounds


class TestBounds(unittest.TestCase):
    
    def testDefaults(self):
        default_lb = 10
        default_ub = 100
        b = bounds.Bounds(default_lb=default_lb,
                          default_ub=default_ub)
        
        keys = ('a', 'askdjn', 'hobos', None, 1, 0.232)
        for t in keys:
            self.assertEquals(default_lb, b.GetLowerBound(t))
            self.assertEquals(default_ub, b.GetUpperBound(t))
        
        lb, ub = b.GetBounds(keys)
        n = len(keys)
        expected_lb = [default_lb] * n
        expected_ub = [default_ub] * n
        self.assertEqual(expected_lb, list(lb))
        self.assertEqual(expected_ub, list(ub))
        
        lb, ub = b.GetLnBounds(keys)
        expected_lb = [np.log(default_lb)] * n
        expected_ub = [np.log(default_ub)] * n
        self.assertEqual(expected_lb, list(lb))
        self.assertEqual(expected_ub, list(ub))
        

    def testSpecifics(self):
        # NOTE(flamholz): values are > 0 on purpose. Logs won't work o.w.
        default_lb = 10
        default_ub = 100
        lb_dict = {'a': 0.0002}
        ub_dict = {'a': 30, 'hobos': 11}
        
        b = bounds.Bounds(lower_bounds=lb_dict,
                          upper_bounds=ub_dict,
                          default_lb=default_lb,
                          default_ub=default_ub)
        
        expected_lb, expected_ub = [], []
        keys = ('a', 'askdjn', 'hobos', None, 1, 0.232)
        for t in keys:
            lb, ub = default_lb, default_ub
            if t in lb_dict:
                lb = lb_dict[t]
            if t in ub_dict:
                ub = ub_dict[t]
            
            expected_ub.append(ub)
            expected_lb.append(lb)
            
            self.assertEquals(lb, b.GetLowerBound(t))
            self.assertEquals(ub, b.GetUpperBound(t))
        
        lb, ub = b.GetBounds(keys)
        self.assertEqual(expected_lb, list(lb))
        self.assertEqual(expected_ub, list(ub))
        
        lb, ub = b.GetLnBounds(keys)
        expected_lb = [np.log(l) for l in expected_lb]
        expected_ub = [np.log(u) for u in expected_ub]
        self.assertEqual(expected_lb, list(lb))
        self.assertEqual(expected_ub, list(ub))
        

def Suite():
    return unittest.makeSuite(TestBounds,'test')
    

if __name__ == '__main__':
    unittest.main()