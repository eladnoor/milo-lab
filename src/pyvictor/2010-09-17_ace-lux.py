import sys, os
from pylab import *
from matplotlib import font_manager
from matplotlib.backends.backend_pdf import PdfPages
from pyvictor.victor_parser import VictorParser

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
if (not os.path.exists('../res/victor')):
    os.mkdir('../res/victor')

name = "elad_100915_ace_lux"
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
t_max = 60

vlegend_k12 = [ ('K12-lux (glu)', 'g', [(0, 0), (0, 2), (0, 4)]), \
            ('K12-lux (succ)', 'b', [(1, 1), (1, 3), (1, 5)]), \
            ('K12-lux (succ+ace)', 'm', [(2, 0), (2, 2), (2, 4)]), \
            ('K12-lux (ace)', 'r', [(3, 1), (3, 3), (3, 5)]) ]

vlegend_ace = [ ('Dace-lux (glu)', 'g', [(0, 6), (0, 8), (0, 10)]), \
            ('Dace-lux (succ)', 'b', [(1, 7), (1, 9), (1, 11)]), \
            ('Dace-lux (succ+ace)', 'm', [(2, 6), (2, 8), (2, 10)]), \
            ('Dace-lux (ace)', 'r', [(3, 7), (3, 9), (3, 11)]) ]

plots.append(('K12', 0, (0, t_max), (3e-2, 1), 'OD', vlegend_k12))
plots.append(('K12', 1, (0, t_max), (1e1, 1e7), 'Lumin', vlegend_k12))
plots.append(('K12', "divide", (0, t_max), (1e2, 1e7), 'Limun/OD', vlegend_k12))

plots.append(('delta-ace', 0, (0, t_max), (3e-2, 1), 'OD', vlegend_ace))
plots.append(('delta-ace', 1, (0, t_max), (1e1, 1e7), 'Lumin', vlegend_ace))
plots.append(('delta-ace', "divide", (0, t_max), (1e2, 1e7), 'Limun/OD', vlegend_ace))

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