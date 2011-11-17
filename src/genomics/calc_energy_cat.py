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
    opt_parser.add_option("-c", "--col",
                          dest="col",
                          help="column to grab")
    return opt_parser

def HandleTag(tag):
    
    if tag in ('photosynthetic',
               'photoautotroph',
               'phototroph'):
        return 'Photoautotroph'
    
    if tag in ('chemolithoautotroph',
               'lithoautotroph'):
        return 'Chemoautotroph'
    
    if tag in ('heterotroph',
               'chemoheterotroph',
               'chemoorganotroph',
               'organotroph'):
        return 'Chemoheterotroph'
    
    if tag in ('methanotroph',
               'methylotroph'):
        return 'C1'
    
    return None


def GetEnergyCategory(energy_tags):
    values = map(HandleTag, energy_tags)
    all_eq = reduce(lambda x,y: x == y, values)
    if all_eq:
        return values[0]
    
    return None
    
    

def Main():
    options, _ = MakeOpts().parse_args(sys.argv)
    assert options.input_filename
    assert options.col

    l = []
    f = open(options.input_filename)
    for row in csv.DictReader(f):
        if options.col in row:
            l.append(row[options.col].lower().split(", "))

    categories = map(GetEnergyCategory, l)
    for cat in categories:
        print cat or ''


if __name__ == '__main__':
    Main()
