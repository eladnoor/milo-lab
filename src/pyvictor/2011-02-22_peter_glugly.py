import sys, os
from pylab import *
from matplotlib import font_manager
from matplotlib.backends.backend_pdf import PdfPages
from pyvictor.victor_parser import VictorParser
from toolbox.util import _mkdir

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

_mkdir('../res/victor')

vp_vec = []
for name in ["peter_glugly"]:
    vp = VictorParser()
    vp.parse_excel("../data/victor/%s.xls" % (name))
    vp_vec.append(vp)

pp = PdfPages('../res/victor/2011-02-22_peter_glugly.pdf')

#rcParams['text.usetex'] = True
rcParams['legend.fontsize'] = 12
#rcParams['font.family'] = 'sans-serif'
#rcParams['font.size'] = 8
#rcParams['lines.linewidth'] = 0.3
#rcParams['lines.markersize'] = 2
#rcParams['figure.figsize'] = [5, 10]
#rcParams['figure.subplot.hspace'] = 0.3
#figure()

plot_growth_rate = False
linewidth = 0.5
fit_window_size = 1.5 # hours
fit_start_threshold = 0.01

plots = [] # (title, victor_index, (t_min, t_max), (y_min, y_max), y_label, 
t_max = 30
OD_min = 0#.046

colors = ['gray', 'red', 'magenta', 'blue', 'cyan', 'green', 'pink', 'orange', 'black']
vlegend = []
vlegend += [('A1-2', 'orange', 'solid', [(0, 0), (0, 1)])]
vlegend += [('A3-4', 'orange', 'dashed', [(0, 2), (0, 3)])]
vlegend += [('B1-2', 'green', 'solid', [(1, 0), (1, 1)])]
vlegend += [('B3-4', 'green', 'dashed', [(1, 2), (1, 3)])]
vlegend += [('B5-6', 'red', 'solid', [(1, 4), (1, 5)])]
vlegend += [('B7-8', 'red', 'dashed', [(1, 6), (1, 7)])]
vlegend += [('B9-10', 'pink', 'solid', [(1, 8), (1, 9)])]
vlegend += [('B11-12', 'pink', 'dashed', [(1, 10), (1, 11)])]
vlegend += [('C1-2', 'blue', 'solid', [(2, 0), (2, 1)])]
vlegend += [('C3-4', 'blue', 'dashed', [(2, 2), (2, 3)])]
vlegend += [('C5-6', 'cyan', 'solid', [(2, 4), (2, 5)])]
vlegend += [('C7-8', 'cyan', 'dashed', [(2, 6), (2, 7)])]
vlegend += [('C5-6', 'black', 'solid', [(2, 8), (2, 9)])]
vlegend += [('C7-8', 'black', 'dashed', [(2, 10), (2, 11)])]
vlegend += [('D1-2', 'magenta', 'solid', [(3, 0), (3, 1)])]
vlegend += [('D3-4', 'magenta', 'dashed', [(3, 2), (3, 3)])]
plots.append(('Rows A-D', (0, t_max), (3e-2, 8e-1), 'OD', vlegend))

vlegend = []
vlegend += [('E1-2', 'green', 'solid', [(4, 0), (4, 1)])]
vlegend += [('E3-4', 'green', 'dashed', [(4, 2), (4, 3)])]
vlegend += [('E5-6', 'red', 'solid', [(4, 4), (4, 5)])]
vlegend += [('E7-8', 'red', 'dashed', [(4, 6), (4, 7)])]
vlegend += [('E5-6', 'orange', 'solid', [(4, 8), (4, 9)])]
vlegend += [('E7-8', 'orange', 'dashed', [(4, 10), (4, 11)])]
vlegend += [('F1-2', 'blue', 'solid', [(5, 0), (5, 1)])]
vlegend += [('F3-4', 'blue', 'dashed', [(5, 2), (5, 3)])]
vlegend += [('F5-6', 'cyan', 'solid', [(5, 4), (5, 5)])]
vlegend += [('F7-8', 'cyan', 'dashed', [(5, 6), (5, 7)])]
vlegend += [('F5-6', 'pink', 'solid', [(5, 8), (5, 9)])]
vlegend += [('F7-8', 'pink', 'dashed', [(5, 10), (5, 11)])]
vlegend += [('G1-2', 'magenta', 'solid', [(6, 0), (6, 1)])]
vlegend += [('G3-4', 'magenta', 'dashed', [(6, 2), (6, 3)])]
plots.append(('Rows E-G', (0, t_max), (3e-2, 8e-1), 'OD', vlegend))

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
            
            try:
                line = plot(time, values, color, linestyle=linestyle, linewidth=linewidth)
            except ValueError:
                raise Exception("the color '%s' is not valid" % color)
                
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