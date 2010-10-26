import sys, os
from pylab import *
from matplotlib import font_manager
from matplotlib.backends.backend_pdf import PdfPages
from pyvictor.victor_parser import VictorParser

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

if (not os.path.exists('../res')):
    os.mkdir('../res')
if (not os.path.exists('../res/victor')):
    os.mkdir('../res/victor')

name = "Elad's OD600_20101025_233"
vp_vec = []
vp = VictorParser()
vp.parse_excel("../data/victor/%s.xls" % (name))
vp_vec.append(vp)

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

fit_window_size = 1.5 # hours
fit_start_threshold = 0.01

plots = [] # (title, victor_index, (t_min, t_max), (y_min, y_max), y_label, 
t_max = 10
OD_min = 0.046

rows_left = ['K12', 'mCherry', '510', '1' ,'2', '8', '10']
rows_right = ['K12', 'prkZ', 'B', 'Csmall', 'D', 'E', 'Z']
colors = ['b', 'm', 'r', 'c', 'k', 'y', 'g']

vlegend = []
for r in xrange(7):
    vlegend += [(rows_left[r] + ' -IPTG', colors[r] + ':', [(r+1, 0), (r+1, 1), (r+1, 2)])]
    vlegend += [(rows_left[r] + ' +IPTG', colors[r], [(r+1, 3), (r+1, 4), (r+1, 5)])]
plots.append(('Left', (0, t_max), (1e-3, 1), 'OD', vlegend))

vlegend = []
for r in xrange(7):
    vlegend += [(rows_right[r] + ' -IPTG', colors[r] + ':', [(r+1, 6), (r+1, 7), (r+1, 8)])]
    vlegend += [(rows_right[r] + ' +IPTG', colors[r], [(r+1, 9), (r+1, 10), (r+1, 11)])]
plots.append(('Right', (0, t_max), (1e-3, 1), 'OD', vlegend))


for (plot_title, t_range, y_range, y_label, data_series) in plots:
    sys.stderr.write("Plotting %s (%s) ... \n" % (plot_title, y_label))
    fig = figure()
    title(plot_title)
    xlabel('Time (hr)')
    ylabel(y_label)
    
    label2legend = {}
    label2line = []
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
            
            line = plot(time, values, color)
            if (label not in label2legend):
                label2line.append((line, label))
                label2legend[label] = label + ", g = "

            try:
                growth_rate = vp.fit_growth(time, values, fit_window_size, fit_start_threshold)
                label2legend[label] += "%.1f  " % growth_rate
            except Exception:
                sys.stderr.write("WARNING: cannot calculate the growth rate in cell (%d, %d)\n" % (row, col))

    rcParams['legend.fontsize'] = 6
    legend([a[0] for a in label2line], [label2legend[a[1]] for a in label2line], loc='lower right')
    yscale('log')
    axis([t_range[0], t_range[1], y_range[0], y_range[1]])
    pp.savefig(fig)

pp.close()