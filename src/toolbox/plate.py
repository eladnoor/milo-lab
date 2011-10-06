#!/usr/bin/python

from toolbox.database import MySQLDatabase
import numpy
import pylab


class TimedMeasurement(object):
    """Measurement with an associated time."""
    
    def __init__(self, time, value):
        self.time = time
        self.value = value


class Plate96(object):
    
    def __init__(self, wells, reading_labels):
        self.wells = wells
        self.reading_labels = set(reading_labels)
    
    @staticmethod
    def FromDatabase(database, exp_id, plate):
        reading_labels = []
        for row in database.Execute("SELECT distinct(reading_label) from tecan_readings WHERE plate='%s' AND exp_id='%s'" % (plate, exp_id)):
            reading_labels.append(str(row[0]))
        
        wells = [[{} for _i in range(12)] for _j in range(8)]
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
        for row in self.wells:
            for well in row:
                measurements = well.get(reading_label)
                readings.append([m.value for m in measurements])
                times.append([m.time for m in measurements])
        
        readings = numpy.vstack(readings)
        times = numpy.vstack(times)
        return times, readings
                
            
if __name__ == '__main__':
    db = MySQLDatabase(host='132.77.80.238', user='ronm', 
                       passwd='a1a1a1', db='tecan')
    exp_id = '2011-10-05 17:55:38'
    plate = '0'
    
    p = Plate96.FromDatabase(db, exp_id, plate)
    times, ods = p.SelectReading('OD600')
    
    pylab.plot(times, ods, 'b.')
    pylab.show()
    
    
    
        