#!/usr/bin/python

"""Calculates the monod curve from a 96 well plate.

Assumes the well labels are substrate concentrations.
"""

import pylab
import sys

from toolbox.database import MySQLDatabase
from toolbox import growth
from toolbox.stats import MeanWithConfidenceInterval
from toolbox.color import ColorMap
from toolbox.plate import Plate96

from optparse import OptionParser


def MakeOpts():
    """Returns an OptionParser object with all the default options."""
    opt_parser = OptionParser()
    opt_parser.add_option("-e", "--experiment_id",
                          dest="experiment_id",
                          help="experiment ID")
    opt_parser.add_option("-p", "--plate_id",
                          dest="plate_id",
                          help="plate ID")
    opt_parser.add_option("-r", "--reading_label",
                          dest="reading_label",
                          default="OD600",
                          help="Reading label")
    opt_parser.add_option("-w", "--window_size",
                          dest="window_size",
                          type='int',
                          default=12,
                          help="Window size for computing the growth rate.")
    opt_parser.add_option("-l", "--lower_bound",
                          dest="lower_bound",
                          type='float',
                          default=0.1,
                          help="Minimum reading to consider valid.")
    opt_parser.add_option("-u", "--upper_bound",
                          dest="upper_bound",
                          type='float',
                          default=0.3,
                          help="Maximum reading to consider valid.")
    return opt_parser


def MaybeFloat(x):
    try:
        return float(x)
    except Exception, e:
        return x
    

def Main():
    options, _ = MakeOpts().parse_args(sys.argv)
    assert options.experiment_id and options.plate_id and options.reading_label
    print 'Reading plate %s from experiment %s' % (options.plate_id, options.experiment_id)
    
    db = MySQLDatabase(host='132.77.80.238', user='ronm', 
                       passwd='a1a1a1', db='tecan')
    p = Plate96.FromDatabase(db, options.experiment_id,
                             options.plate_id)
    times, readings, labels = p.SelectReading(options.reading_label)

    print 'Calculating growth rates'
    lux_calculator = growth.SlidingWindowGrowthCalculator(window_size=options.window_size,
                                                          minimum_level=options.lower_bound,
                                                          maximum_level=options.upper_bound)
    n = labels.size
    colors = ColorMap(set(labels))
    rates = {}
    for i in xrange(n):
        rate, unused_stationary = lux_calculator.CalculateGrowth(times[i,:], readings[i,:])
        scaled_rate = rate * 60 * 60
        label = labels[i]
        rates.setdefault(label, []).append(scaled_rate)
    
    sorted_labels = sorted(rates.keys(), key=MaybeFloat, reverse=False)
    xpts = pylab.arange(len(sorted_labels))
    ticks = xpts + 0.5
    mean_and_err_rates = [MeanWithConfidenceInterval(rates[l])
                          for l in sorted_labels]
    means = [t[0] for t in mean_and_err_rates]
    errs = [t[1] for t in mean_and_err_rates]
    
    for label, mean_rate, err in zip(sorted_labels, means, errs):
        print '%s: %.2g +- %.2g' % (label, mean_rate, err)
    
    pylab.figure()
    pylab.bar(xpts, means)
    pylab.errorbar(ticks, means, yerr=errs, linestyle='None')
    pylab.xticks(ticks, sorted_labels, rotation='vertical')
    pylab.show()
    
    
if __name__ == '__main__':
    Main()