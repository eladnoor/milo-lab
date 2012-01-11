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

def TryFloat(x):
    try:
        return float(x)
    except Exception, e:
        return False


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
    stationaries = {}
    for i in xrange(n):
        label = labels[i]
        if not TryFloat(label):
            continue
        
        rate, stationary = lux_calculator.CalculateGrowth(times[i,:], readings[i,:])
        scaled_rate = rate * 60 * 60
        
        rates.setdefault(label, []).append(scaled_rate)
        stationaries.setdefault(label, []).append(stationary)
    
    print 'Calculating mean growth rate and error'
    f = pylab.figure(0)
    mean_rates = []
    rate_err = []
    for label in rates.iterkeys():
        r = pylab.array(rates[label])
        mean_rate, r_err = MeanWithConfidenceInterval(r)
        mean_rates.append(mean_rate)
        rate_err.append(r_err)
        pct = TryFloat(label)
        if not pct:
            continue
        pct *= 100
        s = pylab.array(stationaries[label]) / pct
        
        s_mean, s_err = MeanWithConfidenceInterval(s)
        if s_mean == 0.0:
            continue
        
        pylab.loglog(pct, s_mean, color=colors[label], marker='.',
                     linestyle='None', figure=f, markersize=20, label=label)
        pylab.errorbar(pct, s_mean, yerr=s_err) 
    pylab.xlabel('Substrate Concentration (%)')
    pylab.ylabel('Yield per Input Concentration (%s/%%)' % (options.reading_label))
    pylab.legend(loc='upper right')
    
    pylab.figure(2)
    concs = pylab.array([float(l) for l in rates.keys()])
    for conc, rate, err in zip(concs, mean_rates, rate_err):
        print '%s: %.2g +- %.2g' % (conc, rate, err)
    
    idx = range(len(concs))
    pylab.semilogx(concs, mean_rates, 'g.')
    pylab.errorbar(concs, mean_rates, yerr=rate_err, fmt=None)
    pylab.xlabel('Substrate Concentration (%)')
    pylab.ylabel('Specific Growth Rate (/hour)')
    pylab.xlim((1.0e-10, concs.max()+1))

    pylab.show()
    
    
if __name__ == '__main__':
    Main()