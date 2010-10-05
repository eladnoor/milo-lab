from victor_parser import VictorParser
from pylab import *
from matplotlib import font_manager
from matplotlib.backends.backend_pdf import PdfPages
import sys, os


def get_data(index, row, col, vp1, vp2):
    (time1, value1) = vp1.get_data(index, row, col)
    (time2, value2) = vp2.get_data(index, row, col)
    
    last_t1 = time1.max()
    time = array([t for t in time1] + [(t + last_t1) for t in time2])
    value = hstack([value1, value2])
    return (time, value)

if (not os.path.exists('../res')):
    os.mkdir('../res')

vp1 = VictorParser()
vp2 = VictorParser()
name = "d100805_elad_LUX"
vp1.parse_excel("../data/%s_1.xls" % name)
vp2.parse_excel("../data/%s_2.xls" % name)
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
t_max = 60

row_labels = ['LUX', 'pZE-lux']
for r in [0, 1]:
    y_label = [('glu 0.1%', 'c:', [(r,0)]), ('glu 0.05%', 'b:', [(r,2)]), ('glu 0.025%', 'k:', [(r,4)]), ('glu 0%', 'r:', [(r,6)]), \
               ('glu 0.1% + 0.1%', 'c', [(r,1)]), ('glu 0.05% + 0.1%', 'b', [(r,3)]), ('glu 0.025% + 0.1%', 'k', [(r,5)]), ('glu 0% + 0.1%', 'r', [(r,7)])]
    plots.append((row_labels[r], 0, (0, t_max), (3e-2, 1), 'OD', y_label))
    plots.append((row_labels[r], 1, (0, t_max), (1e1, 1e7), 'Lumin', y_label))
    plots.append((row_labels[r], "divide", (0, t_max), (1e2, 1e7), 'Limun/OD', y_label))

for (plot_title, victor_index, t_range, y_range, y_label, data_series) in plots:
    sys.stderr.write("Plotting %s (%s) ... \n" % (plot_title, str(victor_index)))
    fig = figure()
    title(plot_title)
    xlabel('Time (hr)')
    ylabel(y_label)
    
    legend_text = []
    for (label, color, cells) in data_series:
        for (row, col) in cells:
            if (victor_index == "divide"):
                (time0, values0) = get_data(0, row, col, vp1, vp2)
                (time1, values1) = get_data(1, row, col, vp1, vp2)
                time = time0
                values = values1/values0                    
            else:
                (time, values) = get_data(victor_index, row, col, vp1, vp2)
            plot(time, values, color)
            #growth_rate = VictorParser.fit_growth(time, values, fit_window_size, fit_start_threshold)
            #legend_text.append("%s (%.0f min)" % (label, log(2)*60/growth_rate))
            legend_text.append("%s" % label)

    rcParams['legend.fontsize'] = 6
    legend(legend_text, loc='upper left')
    yscale('log')
    axis([t_range[0], t_range[1], y_range[0], y_range[1]])
    pp.savefig(fig)

pp.close()