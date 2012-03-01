#!/usr/bin/python

import numpy as np


class FakeThermoData(object):
    
    def GetDGrTagZero_ForModel(self, unused_model):
        return np.matrix([[11.1, -2.1]])
    
    
class FakeInfeasibleThermoData(object):
    
    def GetDGrTagZero_ForModel(self, unused_model):
        return np.matrix([[31.1, 22.1]])