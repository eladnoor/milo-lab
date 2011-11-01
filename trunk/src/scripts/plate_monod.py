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
                          type='int',
                          default=0.1,
                          help="Minimum reading to consider valid.")
    return opt_parser


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
                                                          minimum_level=options.lower_bound)
    n = labels.size
    colors = ColorMap(set(labels))
    rates = {}
    stationaries = {}
    for i in xrange(n):
        rate, stationary = lux_calculator.CalculateGrowth(times[i,:], readings[i,:])
        scaled_rate = rate * 60 * 60
        label = labels[i]
        rates.setdefault(label, []).append(scaled_rate)
        stationaries.setdefault(label, []).append(stationary)
    
    print 'Calculating mean growth rate and error'
    f = pylab.figure(0)
    mean_rates = []
    rate_err = []
    for label in rates.iterkeys():
        r = rates[label]
        mean_rate, err = MeanWithConfidenceInterval(r)
        mean_rates.append(mean_rate)
        rate_err.append(err)
        s = stationaries[label]
        pylab.semilogx(s, r, color=colors[label], marker='.',
                       linestyle='None', figure=f, markersize=20, label=label) 
    pylab.xlabel('Specific Growth Rate (/hour)')
    pylab.ylabel('Stationary Level (%s)' % options.reading_label)
    pylab.legend()
    
    pylab.figure(2)
    concs = pylab.array([float(l) for l in rates.keys()]) * 100
    idx = range(len(concs))
    pylab.semilogx(concs, mean_rates, 'g.')
    pylab.errorbar(concs, mean_rates, yerr=rate_err, fmt=None)
    pylab.xlabel('Substrate Concentration (%)')
    pylab.ylabel('Specific Growth Rate (/hour)')
    pylab.xlim((-0.1, 3.0))

    pylab.show()
    
    
if __name__ == '__main__':
    Main()