#!/usr/bin/python

import numpy as np


class FakeStoichModel(object):
    
    name = 'FakeStoichModel'
    
    def GetStoichiometricMatrix(self):
        return np.matrix([[-1,  0],
                          [ 1, -1],
                          [ 0,  1]])
    
    def GetReactionIDs(self):
        return ['R1', 'R2']
    
    def GetCompoundIDs(self):
        return ['C1', 'C2', 'C3']
    
    def GetFluxes(self):
        return np.array([1.0, 1.0])