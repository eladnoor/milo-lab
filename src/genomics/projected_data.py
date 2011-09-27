#!/usr/bin/python

import logging
import itertools

from genomics.histogram import Histogram


class ProjectedData(object):
    """Project multi-dimensional binary data down to 2d."""
    
    def __init__(self, raw_data):
        self.raw_data = raw_data
    
        self.ind = self.raw_data.ind
        self.ind_vals = self.raw_data.ind_vals
        self.dep_vals = self.ProjectAllBinary(self.raw_data.dep_vals)
    
    def DepCounts(self):
        c = {}
        for dep in self.dep_vals:
            c[dep] = c.get(dep, 0) + 1
        return c
    
    def Iterate(self, filter_values=None):
        ignores = set(filter_values or [])
        for ind, dep in itertools.izip(self.ind_vals, self.dep_vals):
            if dep in ignores:
                continue
            yield ind, dep
        
    def MakeHistogram(self, filter_values=None):
        return Histogram(self, filter_values=filter_values)
        
    @staticmethod
    def ProjectAllBinary(data_dict):
        names = sorted(data_dict.keys())
        lists = [data_dict.get(n) for n in names]
        out = []
        for vals in zip(*lists):
            out.append(ProjectedData.ProjectBinary(vals, names))
        return out
    
    @staticmethod
    def MakeBinary(s):
        return s.lower() == 'true'
    
    @staticmethod
    def ProjectBinary(vals, names):
        if not vals:
            logging.error('Empty vals iterable')
            return None
        
        if len(vals) == 1:
            if ProjectedData.MakeBinary(vals[0]):
                return names[0]
            return 'NONE'
        
        bools = map(ProjectedData.MakeBinary, vals)
        all_true = reduce(lambda x,y: x and y, bools)
        if all_true:
            if len(vals) == 2:
                return 'BOTH'
            return 'ALL'
        
        true_names = [names[i] for i, test in enumerate(bools)
                      if test == True]
        if true_names:
            return ', '.join(true_names)
    
        if len(vals) == 2:
            return 'NEITHER'
        return 'NONE'
