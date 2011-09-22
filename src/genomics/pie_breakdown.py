#!/usr/bin/python

import csv
import itertools
import logging
import pylab
import random
import sys

from genomics.monte_carlo_tester import MonteCarloTester
from genomics.raw_data import RawData
from genomics.row_filterer import RowFilterer
from optparse import OptionParser


def MakeOpts():
    """Returns an OptionParser object with all the default options."""
    opt_parser = OptionParser()
    opt_parser.add_option("-i", "--input_filename",
                          dest="input_filename",
                          help="input CSV")
    opt_parser.add_option("-c", "--col",
                          dest="col",
                          help="Column to break down by")
    opt_parser.add_option("-x", "--filter_cols",
                          dest="filter_cols",
                          help="Columns to filter on, comma separated")
    opt_parser.add_option("-y", "--filter_cols_vals",
                          dest="filter_cols_vals",
                          help="Values to keep, one per column")
    return opt_parser
    

def Main():
    options, _ = MakeOpts().parse_args(sys.argv)
    assert options.input_filename and options.col
    print 'Reading species list from', options.input_filename

    # Set up row filter    
    row_filterer = None
    if options.filter_cols:
        assert options.filter_cols_vals
        filter_cols = map(str.strip, options.filter_cols.split(','))
        filter_cols_vals = map(str.strip, options.filter_cols_vals.split(','))
        row_filterer = RowFilterer(filter_cols, filter_cols_vals)
    
    # Read data and compute histogram.
    data = RawData.FromCsvFile(options.input_filename,
                               options.col,
                               dependent_cols=[],
                               row_filterer=row_filterer)
    ind_counts = data.IndCounts()
    inds = sorted(ind_counts.keys())
    counts = [ind_counts[l] for l in inds]
    labels = [l or 'None given' for l in inds]

    fig1 = pylab.figure(0)
    pylab.pie(counts, labels=labels, autopct='%.1f%%')
    pylab.show()
        



if __name__ == '__main__':
    Main()
