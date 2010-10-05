from victor_parser import VictorParser
from pylab import *
from matplotlib import font_manager
from matplotlib.backends.backend_pdf import PdfPages
import sys, os


def get_data(index, row, col, vp_vec):
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

name = "elad_100903_ace_lux"
vp_vec = []
vp = VictorParser()
vp.parse_excel("../data/%s.xls" % (name))
vp_vec.append(vp)

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
t_max = 35

vlegend = [ ('ace-lux (glu)', 'r', [(0, 0), (0, 2), (0, 4)]), \
            ('K12-lux (glu)', 'r:', [(0, 6), (0, 8), (0, 10)]), \
            ('Dace-lux (succ)', 'b', [(2, 1), (2, 3), (2, 5)]), \
            ('K12-lux (succ)', 'b:', [(2, 7), (2, 9), (2, 11)]), \
            ('Dace-lux (succ+ace)', 'm', [(4, 0), (4, 2), (4, 4)]), \
            ('K12-lux (succ+ace)', 'm:', [(4, 6), (4, 8), (4, 10)]) ] 

plots.append(('M9', 0, (0, t_max), (3e-2, 1), 'OD', vlegend))
plots.append(('M9', 1, (0, t_max), (1e1, 1e7), 'Lumin', vlegend))
plots.append(('M9', "divide", (0, t_max), (1e2, 1e7), 'Limun/OD', vlegend))

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
                (time0, values0) = get_data(0, row, col, vp_vec)
                (time1, values1) = get_data(1, row, col, vp_vec)
                time = time0
                values = values1/values0                    
            else:
                (time, values) = get_data(victor_index, row, col, vp_vec)
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