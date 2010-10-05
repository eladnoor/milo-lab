from victor_parser import VictorParser
from pylab import *
from matplotlib import font_manager
from matplotlib.backends.backend_pdf import PdfPages
import sys, os

if (not os.path.exists('../res')):
    os.mkdir('../res')

vp = VictorParser()
name = "d100805_elad_LUX"
fname = "../data/%s.xls" % name
if (not os.path.exists(fname)):
    raise Exception("Cannot locate the Excel file: " + fname)
vp.parse_excel(fname)
pp = PdfPages('../res/%s.pdf' % name)

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
t_max = 20

row_labels = ['LUX', 'pZE-lux']
for r in [0, 1]:
    y_label = {'glu 0.1\%':('c', [(r,0),(r,1)]), 'glu 0.05\%':('b', [(r,2),(r,3)]), 'glu 0.025\%':('k', [(r,4),(r,5)]), 'glu 0\%':('r', [(r,6),(r,7)])}
    plots.append((row_labels[r], 0, (0, t_max), (3e-2, 7e-1), 'OD', y_label))
    plots.append((row_labels[r], 1, (0, t_max), (1e1, 1e7), 'Lumin', y_label))
    plots.append((row_labels[r], "divide", (0, t_max), (1e2, 1e7), 'Limun/OD', y_label))

for (plot_title, victor_index, t_range, y_range, y_label, data_series) in plots:
    sys.stderr.write("Plotting %s (%s) ... \n" % (plot_title, str(victor_index)))
    fig = figure()
    title(plot_title)
    xlabel('Time (hr)')
    ylabel(y_label)
    
    legend_text = []
    for (label, (color, cells)) in data_series.iteritems():
        for (row, col) in cells:
            if (victor_index == "divide"):
                (time0, values0) = vp.get_data(0, row, col)
                (time1, values1) = vp.get_data(1, row, col)
                time = time0
                values = values1/values0                    
            else:
                (time, values) = vp.get_data(victor_index, row, col)
            plot(time, values, color)
            growth_rate = vp.fit_growth(time, values, fit_window_size, fit_start_threshold)
            legend_text.append("%s (%.0f min)" % (label, log(2)*60/growth_rate))

    legend(legend_text, loc='upper left')
    yscale('log')
    axis([t_range[0], t_range[1], y_range[0], y_range[1]])
    pp.savefig(fig)

pp.close()