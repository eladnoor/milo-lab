#!/usr/bin/python

import logging
import unittest

from toolbox import ambiguous_seq_test
from toolbox import random_seq_test


def main():
    test_modules = (ambiguous_seq_test,
                    random_seq_test)
    
    modules_str = ', '.join(m.__name__ for m in test_modules)
    print 'Running test suites from modules %s' % modules_str
    
    suites = [m.Suite() for m in test_modules]
    alltests = unittest.TestSuite(suites)
    
    runner = unittest.TextTestRunner()
    runner.run(alltests)
    
    
if __name__ == '__main__':
    main() 
    