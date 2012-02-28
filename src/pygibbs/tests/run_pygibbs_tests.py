#!/usr/bin/python

import logging
import unittest

from pygibbs.tests import kegg_compound_test
from pygibbs.tests import kegg_enzyme_test
from pygibbs.tests import pathway_test
from pygibbs.tests import thermo_json_output_test
from pygibbs.tests import group_decomposition_test

from pygibbs.tests.metabolic_modelling import bounds_test
from pygibbs.tests.metabolic_modelling import mtdf_optimizer_test
from pygibbs.tests.metabolic_modelling import protein_optimizer_test
from pygibbs.tests.metabolic_modelling import stoich_model_test
from pygibbs.tests.metabolic_modelling import thermodynamic_data_test


def main():
    test_modules = (kegg_compound_test,
                    kegg_enzyme_test,
                    pathway_test,
                    thermo_json_output_test,
                    group_decomposition_test,
                    bounds_test,
                    mtdf_optimizer_test,
                    protein_optimizer_test,
                    stoich_model_test,
                    thermodynamic_data_test)
    
    modules_str = ', '.join(m.__name__ for m in test_modules)
    logging.info('Running test suites from modules %s' % modules_str)
    
    suites = [m.Suite() for m in test_modules]
    alltests = unittest.TestSuite(suites)
    
    runner = unittest.TextTestRunner()
    runner.run(alltests)
    
    
if __name__ == '__main__':
    main() 
    