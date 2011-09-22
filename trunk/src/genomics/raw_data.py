#!/usr/bin/python

import csv
import random

from genomics.projected_data import ProjectedData


class RawData(object):
    
    def __init__(self, independent_col,
                 independent_vals,
                 dependent_vals,
                 names=None):
        self.ind = independent_col
        self.deps = sorted(dependent_vals.keys())
        self.ind_vals = independent_vals
        self.dep_vals = dependent_vals
        self.names = names
    
    def Project(self):
        return ProjectedData(self)
    
    def IndCounts(self):
        c = {}
        for k in self.ind_vals:
            c[k] = c.get(k, 0) + 1
        return c
    
    @staticmethod
    def FromCsvFile(fname, independent_col, dependent_cols,
                    name_col=None,
                    row_filterer=None):
        ind_vals, dep_vals, names = RawData._ReadData(fname, independent_col,
                                                      dependent_cols,
                                                      name_col=name_col,
                                                      row_filterer=row_filterer)
        return RawData(independent_col, ind_vals, dep_vals,
                       names=names)
    
    @staticmethod
    def _ReadData(fname, independent_col, dependent_cols,
                  name_col=None,
                  row_filterer=None):
        f = open(fname)
        r = csv.DictReader(f)
        ind = []
        d = dict((c,[]) for c in dependent_cols)
        names = []

        for row in r:
            if row_filterer and not row_filterer.Keep(row):
                continue

            if name_col:
                names.append(row.get(name_col, None))
                
            ind.append(row.get(independent_col, None))
            store_val = lambda c: d[c].append(row.get(c, None))
            map(store_val, dependent_cols)

        f.close()
        return ind, d, names or None
    
    def Shuffle(self):
        new_ind = list(self.ind_vals)
        
        new_dep = dict((key, list(value))
                       for key, value in self.dep_vals.iteritems())
        for key in new_dep.keys():
            l = new_dep[key]
            random.shuffle(l)
        
        return RawData(self.ind, new_ind, new_dep)
    