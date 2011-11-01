#!/usr/bin/python

from toolbox.database import MySQLDatabase
from toolbox import growth
from toolbox.stats import MeanWithConfidenceInterval
from toolbox.color import ColorMap
import numpy
import pylab
import sys

from optparse import OptionParser


class TimedMeasurement(object):
    """Measurement with an associated time."""
    
    def __init__(self, time, value):
        self.time = time
        self.value = value


class Well(dict):
    """Class containing per-well data."""
        
    def __init__(self, row=None, col=None, label=None):
        self.row = row
        self.col = col
        self.label = label
    
    def update(self, row, col, label):
        self.row = row
        self.col = col
        self.label = label


class Plate96(object):
    
    def __init__(self, wells, reading_labels):
        self.wells = wells
        self.reading_labels = set(reading_labels)
    
    @staticmethod
    def FromDatabase(database, exp_id, plate):
        reading_labels = []
        for row in database.Execute("SELECT distinct(reading_label) from tecan_readings WHERE plate='%s' AND exp_id='%s'" % (plate, exp_id)):
            reading_labels.append(str(row[0]))
        
        wells = [[Well() for _i in range(12)] for _j in range(8)]
        sql = "SELECT row, col, label FROM tecan_labels WHERE plate='%s' AND exp_id='%s'" % (plate, exp_id)
        well_labels = database.Execute(sql)
        for row, col, label in well_labels:
            wells[row][col].update(row, col, label)
        
        readings = database.Execute("SELECT row, col, reading_label, time, measurement from tecan_readings WHERE plate='%s' AND exp_id='%s'" % (plate, exp_id))
        for r in readings:
            row, col, reading_label, time, measurement = r            
            
            m = TimedMeasurement(time, measurement)
            wells[row][col].setdefault(reading_label, []).append(m)
        
        return Plate96(wells, reading_labels)
    
    def SelectReading(self, reading_label):
        if reading_label not in self.reading_labels:
            return []
        
        readings = []
        times = []
        labels = []
        for row in self.wells:
            for well in row:
                labels.append(well.label)
                measurements = well.get(reading_label)
                readings.append([m.value for m in measurements])
                times.append([m.time for m in measurements])
                
        
        readings = numpy.vstack(readings)
        times = numpy.vstack(times)
        return times, readings, numpy.array(labels)
                

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
    pylab.loglog(concs, mean_rates, 'g.')
    pylab.errorbar(concs, mean_rates, yerr=rate_err, fmt=None)
    pylab.xlabel('Substrate Concentration (%)')
    pylab.ylabel('Specific Growth Rate (/hour)')
    pylab.xlim((-0.1, 3.0))

    pylab.show()
    
    
if __name__ == '__main__':
    Main()
    
    
    
    
        