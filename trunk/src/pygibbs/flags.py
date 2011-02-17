#!/usr/bin/python

from pygibbs.thermodynamic_constants import R, default_I, default_pH
from pygibbs.thermodynamic_constants import default_pMg, default_T
from optparse import OptionParser


def MakeOpts():
    """Returns an OptionParser object with all the default options."""
    opt_parser = OptionParser()
    opt_parser.add_option("-c", "--median_concentration", dest="c_mid",
                          default=1e-3, type="float",
                          help="Median concentation, Molar")
    opt_parser.add_option("-p", "--ph", dest="ph", type="float",
                          default=default_pH, help="The pH")
    opt_parser.add_option("-m", "--pmg", dest="pmg", type="float",
                          default=default_pMg, help="The pMg")
    opt_parser.add_option("-i", "--ionic_strength", dest="i_s", type="float",
                          default=default_pMg, help="The ionic strength, mM")
    opt_parser.add_option("-g", "--ignore_cofactors", action="store_true",
                          dest="ignore_cofactors", default=False,
                          help="If True, don't fix co-factor concentrations.")
    
    return opt_parser