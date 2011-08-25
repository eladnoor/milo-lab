#!/usr/bin/python

import StringIO
import unittest

from genomics import row_filterer


class TestRowFilterer(unittest.TestCase):
    
    def testSimpleFilter(self):
        filterer = row_filterer.RowFilterer(['a'], [1])
        test_rows = (({'a':1}, True),
                     ({'a':2}, False),
                     ({'b':1}, False),
                     ({'b':2}, False),
                     ({'a':1, 'b':2}, True),
                     ({'a':2, 'b':1}, False),
                     ({'a':2, 'b':2}, False))
        
        for row, keep_val in test_rows:
            self.assertEquals(keep_val, filterer.Keep(row))
    
    def testTwoColumnFilter(self):
        filterer = row_filterer.RowFilterer(['a', 'b'], [1, 2])
        test_rows = (({'a':1}, True),
                     ({'a':2}, False),
                     ({'b':1}, False),
                     ({'b':2}, True),
                     ({'a':1, 'b':2}, True),
                     ({'a':2, 'b':1}, False),
                     ({'a':2, 'b':2}, True))
        
        for row, keep_val in test_rows:
            self.assertEquals(keep_val, filterer.Keep(row))
            

def Suite():
    suites = (unittest.makeSuite(TestRowFilterer,'test'),)
    return unittest.TestSuite(suites)


if __name__ == '__main__':
    unittest.main()