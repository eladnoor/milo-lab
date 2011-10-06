#!/usr/bin/python

from toolbox.database import MySQLDatabase
from toolbox import growth
import numpy
import pylab


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
                
            
if __name__ == '__main__':
    db = MySQLDatabase(host='132.77.80.238', user='ronm', 
                       passwd='a1a1a1', db='tecan')
    exp_id = '2011-10-05 17:55:38'
    plate = '2'
    
    p = Plate96.FromDatabase(db, exp_id, plate)
    times, ods, labels = p.SelectReading('OD600')
    blanks = pylab.find(labels == 'BLANK')
    promoterless = pylab.find(labels == 'AZ01 F3 U66')


    f1 = pylab.figure(0)
    pylab.plot(times[blanks, :], ods[blanks, :], 'b.', figure=f1)
    pylab.plot(times[promoterless, :], ods[promoterless, :], 'g.', figure=f1)    
    

    f2 = pylab.figure(1)
    pylab.plot(times, numpy.log(ods), 'b.', figure=f2)
    
    
    
    calculator = growth.SlidingWindowGrowthCalculator()
    rows, _ = times.shape
    rates = []
    stationaries = []
    for i in xrange(rows):
        rate, stationary = calculator.CalculateGrowth(times[i, :], ods[i, :])
        rates.append(rate)
        stationaries.append(stationary)
    
    f3 = pylab.figure(2)
    pylab.hist(rates, figure=f3)
    
    f4 = pylab.figure(3)
    pylab.hist(stationaries, figure=f3)
    pylab.show()
    
    
    
    
    
        