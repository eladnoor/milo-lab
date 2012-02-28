#!/usr/bin/python

import numpy as np


class FakeThermoData(object):
    
    def GetDGrTagZero_ForModel(self, unused_model):
        return np.array([11.1, -2.1])