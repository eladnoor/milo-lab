#!/usr/bin/python

"""Calculates the monod curve from a 96 well plate.

Assumes the well labels are substrate concentrations.
"""

import pylab
import sys

from toolbox.database import MySQLDatabase
from toolbox import growth
from toolbox import smoothing
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
    opt_parser.add_option("-p", "--plate_ids",
                          dest="plate_ids",
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


def MeanWithConfidenceIntervalDict(d):
    return dict((k, MeanWithConfidenceInterval(v))
                for k,v in d.iteritems())


def Main():
    options, _ = MakeOpts().parse_args(sys.argv)
    assert options.experiment_id and options.plate_ids and options.reading_label
    plates = map(str.strip, options.plate_ids.split(','))
    
    print 'Reading plates %s from experiment %s' % (', '.join(plates),
                                                    options.experiment_id)
    
    db = MySQLDatabase(host='hldbv02', user='ronm', 
                       passwd='a1a1a1', db='tecan')

    print 'Calculating growth rates'
    growth_calc = growth.SlidingWindowGrowthCalculator(window_size=options.window_size,
                                                       minimum_level=options.lower_bound,
                                                       maximum_level=options.upper_bound)
    
    plate_names = {'1': 'Glucose',
                   '2': 'Gluconate'}
    colormap = {'1': 'k',
                '2': 'r'}
    f1 = pylab.figure(0)
    for plate_id in plates:
        p = Plate96.FromDatabase(db, options.experiment_id, plate_id)
        rates, unused_stationaries = growth_calc.CalculatePlateGrowth(
            p, options.reading_label)
        mean_rates = MeanWithConfidenceIntervalDict(rates)
        
        means = []
        errors = []
        concs = []
        for label in mean_rates.keys():
            conc = TryFloat(label)
            if conc is False:
                continue
            
            mean, error = mean_rates[label]
            means.append(mean)
            errors.append(error)
            concs.append(conc)
        
        means = pylab.array(means)
        errors = pylab.array(errors)
        concs = pylab.array(concs)
        
        max_mean = max(means)
        norm_means = means / max_mean
        norm_errors = errors / max_mean
        
        label = plate_names[plate_id]
        color = colormap[plate_id]
        """
        pylab.subplot(121)
        pylab.plot(concs, norm_means, color=color, linestyle='None',
                   marker='.', label=label)
        pylab.errorbar(concs, norm_means, yerr=norm_errors, ecolor=color,
                       fmt=None)
        """
        
        #pylab.subplot(122)
        pcts = concs*100
        idx = pylab.find(pcts > 1e-4)
        my_pcts = pcts[idx]
        my_means = means[idx]
        my_errs = errors[idx]
        
        order = pylab.argsort(my_pcts)
        my_pcts = my_pcts[order]
        my_means = my_means[order]
        my_errs = my_errs[order]
        
        pylab.plot(my_pcts, my_means, color=color, linestyle='None',
                   linewidth=4, marker='.', markersize=15, label=label)
        pylab.errorbar(my_pcts, my_means, yerr=my_errs, ecolor=color,
                       fmt=None, linewidth=1)
        
        smoothed = smoothing.WeightedAverageSmoother(pylab.log(my_pcts), my_means,
                                                     sigma=0.7)
        log_xs = pylab.arange(pylab.log(1e-4), pylab.log(2.2), 1e-3)
        xs = pylab.exp(log_xs)
        ys = smoothed(log_xs)
        pylab.plot(xs, ys, color=color, linewidth=3, linestyle='--')
        
        
    
    """            
    pylab.subplot(121)
    pylab.xlabel('CAP Concentration (fraction of standard concentration)')
    pylab.ylabel('Relative Specific Growth Rate (/hour)')
    pylab.xlim(-0.1,0.2)
    """
    
    """
    pylab.subplot(122)
    pylab.xlabel('CAP Concentration (fraction of standard concentration)')
    pylab.ylabel('Absolute Specific Growth Rate (/hour)')
    """
    
    pylab.xscale('log')
    pylab.xlabel('Substrate concentration (m/v %)')
    pylab.ylabel('Specific Growth Rate (/hour)')
    
    #pylab.xlim(-0.1,0.2)
    pylab.legend(loc='upper left')
    pylab.show()
    
    
if __name__ == '__main__':
    Main()