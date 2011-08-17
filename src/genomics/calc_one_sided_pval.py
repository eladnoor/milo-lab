#!/usr/bin/python

import csv
import itertools
import logging
import sys

from collections import Counter
from optparse import OptionParser
from colormap import ColorMap
from numpy import random


def MakeOpts():
    """Returns an OptionParser object with all the default options."""
    opt_parser = OptionParser()
    opt_parser.add_option("-i", "--input_filename",
                          dest="input_filename",
                          help="input CSV")
    opt_parser.add_option("-r", "--row",
                          dest="row",
                          help="Row name")
    opt_parser.add_option("-t", "--num_trials",
                          dest="num_trials",
                          type="int",
                          default=10000,
                          help="Number of trials")
    return opt_parser


def Main():
    options, _ = MakeOpts().parse_args(sys.argv)
    assert options.input_filename
    assert options.row
    print 'Reading data', options.input_filename
    
    f = open(options.input_filename)
    r = csv.DictReader(f)
    fieldnames = list(r.fieldnames)
    all_col_keys = set(fieldnames[1:])
    
    data_dict = {}
    col_key_total = 0
    all_row_keys = set()
    for row in r:
        key = row.pop('', None)
        if not key:
            continue
        
        all_row_keys.add(key)
        int_dict = dict((k, int(v)) for k, v in row.iteritems())
        data_dict[key] = int_dict
    
    row_key_total = sum(data_dict[options.row].values())
    
    print '%s has %d entries total' % (options.row, row_key_total)
    
    trials = options.num_trials
    print 'Sampling %d times' % trials
    
    total_ge = dict((key, 0) for key in all_col_keys)
    total_le = dict((key, 0) for key in all_col_keys)
    for _ in xrange(trials):
        population = list(all_col_keys)
        samples = random.randint(0, len(population), row_key_total)
        counts = Counter(samples)
        sampled_counts = dict((population[i], count)
                              for i, count in counts.iteritems())
        
        
        for col_key in all_col_keys:
            my_val = data_dict[options.row][col_key]
            sample_val = sampled_counts[col_key]
            
            if sample_val >= my_val:
                total_ge[col_key] = total_ge.get(col_key, 0) + 1
            if sample_val <= my_val:
                total_le[col_key] = total_le.get(col_key, 0) + 1
                    
            

    for col_key in all_col_keys:
        print 'Observed count for (%s, %s): %d' % (options.row,
                                                   col_key,
                                                   data_dict[options.row][col_key])
        
        print 'Samples lower:', total_le[col_key]
        print 'Samples higher:', total_ge[col_key]
        print 'Lower p-value: %.2g' % (float(total_le[col_key]) / float(trials))
        print 'Higher p-value: %.2g' % (float(total_ge[col_key]) / float(trials))
    
    
if __name__ == '__main__':
    Main()