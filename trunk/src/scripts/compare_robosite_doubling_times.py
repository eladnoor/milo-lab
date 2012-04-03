#!/usr/bin/python

import csv
import pylab

from toolbox import stats


def ReadAndAvgDoublingTimes(filename):
    d = {}
    f = open(filename)
    for row in csv.reader(f):
        key = row[0]
        doubling_time = float(row[-1])
        d.setdefault(key, []).append(doubling_time)
    
    out = {}
    for key, dts in d.iteritems():
        out[key] = stats.MeanWithConfidenceInterval(dts)
        
    return out


dts1 = ReadAndAvgDoublingTimes('/home/flamholz/Desktop/DoublingEst_WT.txt')
dts2 = ReadAndAvgDoublingTimes('/home/flamholz/Desktop/DoublingEst_Adapted.txt')


keys = sorted(dts1.keys())
vals1 = [dts1[k][0] for k in keys]
vals2 = [dts2[k][0] for k in keys]
errs1 = [dts1[k][1] for k in keys]
errs2 = [dts2[k][1] for k in keys]

pylab.figure()
pylab.xlabel('WT Doubling Time')
pylab.ylabel('Adapted Doubling Time')
pylab.errorbar(vals1, vals2, yerr=errs2, xerr=errs1, fmt='g.', ecolor='b')
#pylab.errorbar(x, y, yerr, xerr, fmt, ecolor, elinewidth, capsize, barsabove, lolims, uplims, xlolims, xuplims, hold)
for i, k in enumerate(keys):
    pylab.text(vals1[i], vals2[i], k)

pylab.show()

    

