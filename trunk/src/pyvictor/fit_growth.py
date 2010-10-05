"""
    time       - a vector of the time of each measurement (in hours)
    cell_count - a vector of the cell-count at each time point
    window     - the size of the window for which to calculate the growth rates (in hours)
    start_time - the time point after which the input data can be trusted

    Note: the cell_count can be any measure that is proportional to the
    number of cells. For example, Optical Density (in OD), or level of
    luminescence (in CPS).
"""
from victor_parser import VictorParser
import pylab

if (__name__ == "__main__"):
    vp = VictorParser()
    fname = "../data/Elad's OD600_20100701_188.xls"
    vp.parse_excel(fname)
    (time, cell_count) = vp.get_data(0, 3, 0)

    pylab.figure()
    pylab.subplot(1,2,1)
    print "Logistic fit: gr = %.2f" % vp.fit_growth2(time, cell_count, plot_figure=True)
    pylab.subplot(1,2,2)
    print "linear fit in logscale: gr = %.2f" %  vp.fit_growth(time, cell_count, window_size=1.5, start_threshold=0.01, plot_figure=True)
    pylab.show()