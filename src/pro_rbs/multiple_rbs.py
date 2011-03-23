#!/usr/bin/python

import csv
import logging
import os
import pylab
import sys

from optparse import OptionParser
from pro_rbs import rbs_util
from toolbox import ambiguous_seq


# NOTE(flamholz): set this to the location of your NUPACK install.
os.environ['NUPACKHOME'] = '/home/flamholz/Downloads/nupack3.0/'


def MakeOpts():
    """Returns an OptionParser object with all the default options."""
    opt_parser = OptionParser()
    opt_parser.add_option("-i", "--rbs_file", dest="rbs_filename",
                          default="pro_rbs/rbs_v1.csv",
                          help=("The name of the file containing all the RBS sequences"
                                " including spacers, start codon, and (optionally) his_tag."))
    return opt_parser


def main():
    options, _ = MakeOpts().parse_args(sys.argv)
    print 'RBS filename:', options.rbs_filename
    
    f = open(options.rbs_filename)
    dg_by_seq = {}
    failure_count = 0
    
    for row in csv.reader(f):
        rbs = row[0].upper()
        print 'Testing RBS', rbs
        try:
            dG = rbs_util.TryCalcDG(rbs)
            print 'Predicted dG %.2g' % rbs_util.TokJ(dG)
            print 'Predicted K_binding %f' % rbs_util.PredictedRate(dG)
            dg_by_seq[rbs] = dG
        except Exception, e:
            logging.error(e)
            failure_count += 1
            continue
    
    print 'Tested %d RBS options.' % len(dg_by_seq)
    print '%d RBS options failed.' % failure_count
    
    binding_dGs = dg_by_seq.values()
    log_rates = map(rbs_util.LogRate, binding_dGs)
    pylab.hist(log_rates)
    pylab.xlim((-4, 4))
    pylab.xlabel('log10(K_binding)')
    pylab.ylabel('Number of RBS')
    pylab.show()
    

if __name__ == '__main__':
    main()
    
