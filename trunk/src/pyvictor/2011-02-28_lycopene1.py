import sys, os
from pylab import *
from matplotlib import font_manager
from matplotlib.backends.backend_pdf import PdfPages
from pyvictor.victor_parser import VictorParser
from toolbox.util import _mkdir
from pyvictor.util import RowCol2String

def get_data(index, row, col, vp_vec):
    """
        When the experimental data is broken into more than one XLS sheet, this method
        concatenates the data into one series and returns it as if it was from one source.
    """
    time_array = array([])
    value_array = array([])
    last_t = 0
    for vp in vp_vec:
        times, values = vp.get_data(index, row, col)
        time_array = hstack([time_array, times + last_t])
        value_array = hstack([value_array, values])
        if len(time_array) > 0:
            last_t = time_array.max()

    return time_array, value_array

_mkdir('../res/victor')

vp_vec = []
for name in ["OD600 20110302_lycopene1"]:
    vp = VictorParser()
    vp.parse_excel("../data/victor/%s.xls" % (name))
    vp_vec.append(vp)

pp = PdfPages('../res/victor/2011-02-28_lycopene1.pdf')

#rcParams['text.usetex'] = True
rcParams['legend.fontsize'] = 12
#rcParams['font.family'] = 'sans-serif'
#rcParams['font.size'] = 8
#rcParams['lines.linewidth'] = 0.3
#rcParams['lines.markersize'] = 2
#rcParams['figure.figsize'] = [5, 10]
#rcParams['figure.subplot.hspace'] = 0.3
#figure()

plot_growth_rate = True
linewidth = 1
fit_window_size = 1.5 # hours
fit_start_threshold = 0.01

plots = [] # (title, victor_index, (t_min, t_max), (y_min, y_max), y_label, 
t_max = 40
OD_min = 0.046

colors = ['green', 'cyan', 'blue', 'red', 'magenta', 'orange', 'pink', 'orange', 'black']
vlegend = []
for i in xrange(2):
    for r in xrange(5):
        vlegend = []
        vlegend += [(RowCol2String(r, 5*i+j), colors[j], 'solid', [(r, 5*i+j)]) for j in xrange(3)]
        plots.append((RowCol2String(r, 5*i) + '-' + RowCol2String(r, 5*i+2), (0, t_max), (1e-3, 0.4), 'OD', vlegend))

for (plot_title, t_range, y_range, y_label, data_series) in plots:
    sys.stderr.write("Plotting %s (%s) ... \n" % (plot_title, y_label))
    fig = figure()
    title(plot_title)
    xlabel('Time (hr)')
    ylabel(y_label)
    
    label2legend = {}
    label2line = []
    for (label, color, linestyle, cells) in data_series:
        for (row, col) in cells:
            if (y_label == 'Limun/OD'):
                (time0, values0) = get_data(0, row, col, vp_vec)
                (time1, values1) = get_data(1, row, col, vp_vec)
                time = time0
                values = values1/(values0 - OD_min)                    
            elif (y_label == 'OD'):
                (time, values) = get_data(0, row, col, vp_vec)
                values -= OD_min
            elif (y_label == 'Lumin'):
                (time, values) = get_data(1, row, col, vp_vec)
            else:
                raise Exception("unrecognised Y label: " + y_label)
            
            line = plot(time, values, color, linestyle=linestyle, linewidth=linewidth)
            if (label not in label2legend):
                label2line.append((line, label))
                label2legend[label] = label
                if plot_growth_rate:
                    label2legend[label] += ", T(min) = "
            
            if plot_growth_rate:
                try:
                    growth_rate = vp.fit_growth(time, values, fit_window_size, fit_start_threshold)
                except Exception:
                    sys.stderr.write("WARNING: cannot calculate the growth rate in cell (%d, %d)\n" % (row, col))
                if (growth_rate > 1e-10):
                    label2legend[label] += "%.0f  " % (60.0 * log(2.0) / growth_rate)
                else:
                    label2legend[label] += "0  "

    rcParams['legend.fontsize'] = 6
    legend([a[0] for a in label2line], [label2legend[a[1]] for a in label2line], loc='lower right')
    yscale('log')
    axis([t_range[0], t_range[1], y_range[0], y_range[1]])
    pp.savefig(fig)

pp.close()