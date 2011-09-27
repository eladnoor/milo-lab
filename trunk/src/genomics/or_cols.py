#!/usr/bin/python

import csv
import sys

from optparse import OptionParser


def MakeOpts():
    """Returns an OptionParser object with all the default options."""
    opt_parser = OptionParser()
    opt_parser.add_option("-i", "--input_filename",
                          dest="input_filename",
                          help="input CSV")
    opt_parser.add_option("-c", "--cols",
                          dest="cols",
                          help="comma-separated column names")
    opt_parser.add_option("-a", "--and_mode",
                          dest="and_mode",
                          default=False,
                          action="store_true",
                          help="turn on AND mode")
    
    return opt_parser

def ToBool(string):
    if not string:
        return False
    
    l = string.lower()
    return l == 'true'

def Main():
    options, _ = MakeOpts().parse_args(sys.argv)
    assert options.input_filename and options.cols

    colnames = map(str.strip, options.cols.split(','))
    

    r = csv.DictReader(open(options.input_filename))
    for row in r:
        c = [ToBool(row[c]) for c in colnames]
        bool = reduce(lambda x, y: x or y, c)
        if options.and_mode:
            bool = reduce(lambda x, y: x and y, c)
        print bool
        


if __name__ == '__main__':
    Main()