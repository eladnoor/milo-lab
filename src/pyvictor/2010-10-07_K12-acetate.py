import sys, os
from pylab import *
from matplotlib import font_manager
from matplotlib.backends.backend_pdf import PdfPages
from pyvictor.victor_parser import VictorParser
from toolbox import util

def get_data(index, row, col, vp_vec):
    """
        When the experimental data is broken into more than one XLS sheet, this method
        concatenates the data into one series and returns it as if it was from one source.
    """
    time_array = array([])
    value_array = array([])
    last_t = 0
    for vp in vp_vec:
        (time, value) = vp.get_data(index, row, col)
        time_array = hstack([time_array, time + last_t])
        value_array = hstack([value_array, value])
        if (len(time_array) > 0):
            last_t = time_array.max()

    return (time_array, value_array)

name = "2010-10-07_K12-acetate"
vp_vec = []
vp = VictorParser()
vp.parse_excel("../data/victor/%s.xls" % (name))
vp_vec.append(vp)

util._mkdir('../res/victor')
pp = PdfPages('../res/victor/%s.pdf' % name)

#rcParams['text.usetex'] = True
rcParams['legend.fontsize'] = 12
#rcParams['font.family'] = 'sans-serif'
#rcParams['font.size'] = 8
#rcParams['lines.linewidth'] = 0.3
#rcParams['lines.markersize'] = 2
#rcParams['figure.figsize'] = [5, 10]
#rcParams['figure.subplot.hspace'] = 0.3
#figure()

fit_window_size = 5 # hours
fit_start_threshold = 0.002

plots = [] # (title, victor_index, (t_min, t_max), (y_min, y_max), y_label, 
t_max = 30
OD_min = 0.046

rows = ['glucose', 'ribose', 'glycolate', 'glycine', 'succinate', 'glycerol', 'mannose', 'pyruvate']
colors = ['g', 'b', 'm', 'r', 'c', 'k', 'b:', 'g:']

vlegend = []
for r in xrange(8):
    vlegend += [(rows[r], colors[r], [(r, 0), (r, 1)])]
plots.append(('0.05% sugar', (0, t_max), (1e-3, 1), 'OD', vlegend))

vlegend = []
for r in xrange(8):
    vlegend += [(rows[r], colors[r], [(r, 2), (r, 3)])]
plots.append(('0.025% sugar', (0, t_max), (1e-3, 1), 'OD', vlegend))

vlegend = []
for r in xrange(8):
    vlegend += [(rows[r], colors[r], [(r, 4), (r, 5)])]
vlegend += [('nothing', 'y', [(1, 7), (1, 8), (2, 7), (2, 8)])]
vlegend += [('acetate', 'm:', [(5, 7), (5, 8), (6, 7), (6, 8)])]
plots.append(('0.025% sugar + 0.025% acetate', (0, t_max), (1e-3, 1), 'OD', vlegend))


for (plot_title, t_range, y_range, y_label, data_series) in plots:
    sys.stderr.write("Plotting %s (%s) ... \n" % (plot_title, y_label))
    fig = figure()
    title(plot_title)
    xlabel('Time (hr)')
    ylabel(y_label)
    
    label2line = {}
    for (label, color, cells) in data_series:
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
            #growth_rate = vp.fit_growth(time, values, fit_window_size, fit_start_threshold)
            #label += " (%.2f 1/h)" % growth_rate
            label2line[label] = plot(time, values, color)

    rcParams['legend.fontsize'] = 6
    legend(label2line.values(), label2line.keys(), loc='upper left')
    yscale('log')
    axis([t_range[0], t_range[1], y_range[0], y_range[1]])
    pp.savefig(fig)

pp.close()