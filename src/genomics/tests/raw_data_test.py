#!/usr/bin/python

import StringIO
import unittest

from genomics.raw_data import RawData


class TestRawData(unittest.TestCase):
    
    def testEmptyData(self):
        # Should be fine
        raw_data = RawData(None,[],{})
        raw_data.Shuffle()
        raw_data.Project()
    
    def MakeSimpleRawData(self):
        ind_key = 'a'
        ind_vals = [1,2,3,4]
        dep_vals = {'b': ['true', 'false', 'true', 'false'],
                    'c': ['true', 'true', 'true', 'false']}
        return RawData(ind_key, ind_vals, dep_vals)
        
    
    def testSimpleDataShuffle(self):
        raw_data = self.MakeSimpleRawData()
        shuffled = raw_data.Shuffle()
        
        self.assertEquals(raw_data.ind, shuffled.ind)
        self.assertEquals(raw_data.ind_vals, shuffled.ind_vals)
        self.assertEquals(set(raw_data.deps), set(shuffled.deps))
        
        for key in raw_data.deps:
            self.assertEquals(set(raw_data.dep_vals[key]),
                              set(shuffled.dep_vals[key]))
            
    def testSimpleDataProject(self):
        raw_data = self.MakeSimpleRawData()        
        projected = raw_data.Project()
        
        ind_vals = set()
        for ind, dep in projected.Iterate():
            ind_vals.add(ind)
        
        self.assertEquals(set(raw_data.ind_vals), ind_vals)
            

def Suite():
    suites = (unittest.makeSuite(TestRawData,'test'),)
    return unittest.TestSuite(suites)


if __name__ == '__main__':
    unittest.main()