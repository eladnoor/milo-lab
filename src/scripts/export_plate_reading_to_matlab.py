#!/usr/bin/python

"""Calculates the monod curve from a 96 well plate.

Assumes the well labels are substrate concentrations.
"""

import pylab
import sys

import scipy.io as sio
import numpy as np

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
    opt_parser.add_option("-o", "--matfile_name",
                          dest="matfile_name",
                          help="Name of matlab file to create.")
    return opt_parser


def Main():
    options, _ = MakeOpts().parse_args(sys.argv)
    assert options.experiment_id and options.plate_id and options.reading_label
    assert options.matfile_name
    print 'Reading plate %s from experiment %s' % (options.plate_id, options.experiment_id)
    
    db = MySQLDatabase(host='132.77.80.238', user='ronm', 
                       passwd='a1a1a1', db='tecan')
    p = Plate96.FromDatabase(db, options.experiment_id,
                             options.plate_id)
    times, readings, labels = p.SelectReading(options.reading_label)
    
    labels = pylab.array(labels, dtype=np.object)
    
    print 'Writing data to %s' % options.matfile_name
    out_data = {'times': times,
                'labels': labels,
                options.reading_label: readings}
    sio.savemat(options.matfile_name, out_data, oned_as='column')
    
if __name__ == '__main__':
    Main()