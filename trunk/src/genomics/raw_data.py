#!/usr/bin/python

import csv
import random

from genomics.projected_data import ProjectedData


class RawData(object):
    
    def __init__(self, independent_col,
                independent_vals,
                dependent_vals):
        self.ind = independent_col
        self.deps = sorted(dependent_vals.keys())
        self.ind_vals = independent_vals
        self.dep_vals = dependent_vals
    
    def Project(self):
        return ProjectedData(self)
    
    @staticmethod
    def FromCsvFile(fname, independent_col, dependent_cols, row_filterer=None):
        ind_vals, dep_vals = RawData._ReadData(fname, independent_col,
                                               dependent_cols,
                                               row_filterer=row_filterer)
        return RawData(independent_col, ind_vals, dep_vals)
    
    @staticmethod
    def _ReadData(fname, independent_col, dependent_cols, row_filterer=None):
        f = open(fname)
        r = csv.DictReader(f)
        ind = []
        d = dict((c,[]) for c in dependent_cols)

        for row in r:
            if row_filterer and not row_filterer.Keep(row):
                continue
            
            ind.append(row.get(independent_col, None))
            
            store_val = lambda c: d[c].append(row.get(c, None))
            map(store_val, dependent_cols)

        f.close()
        return ind, d
    
    def Shuffle(self):
        new_ind = list(self.ind_vals)
        
        new_dep = dict((key, list(value))
                       for key, value in self.dep_vals.iteritems())
        for key in new_dep.keys():
            l = new_dep[key]
            random.shuffle(l)
        
        return RawData(self.ind, new_ind, new_dep)
    