#!/usr/bin/python

import logging
import unittest

from pygibbs.tests import kegg_compound_test
from pygibbs.tests import kegg_enzyme_test
from pygibbs.tests import thermo_json_output_test


def main():
    test_modules = (kegg_compound_test,
                    kegg_enzyme_test,
                    thermo_json_output_test)
    
    modules_str = ', '.join(m.__name__ for m in test_modules)
    print 'Running test suites from modules %s' % modules_str
    
    suites = [m.Suite() for m in test_modules]
    alltests = unittest.TestSuite(suites)
    
    runner = unittest.TextTestRunner()
    runner.run(alltests)
    
    
if __name__ == '__main__':
    main() 
    