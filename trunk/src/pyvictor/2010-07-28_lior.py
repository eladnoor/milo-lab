from victor_parser import VictorParser
from pylab import *
from matplotlib import font_manager
from matplotlib.backends.backend_pdf import PdfPages
import sys, os

if (not os.path.exists('../res')):
    os.mkdir('../res')

vp = VictorParser()
name = "RBS-GFP 28-07-10"
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

plot_types = [(0, "OD", 0.001, 1), (1, "GFP", 1e3, 1e6), ("divide", "GFP/OD", 1e4, 2e6)]

columns = [('A', [8], 25),('B', [9], 25),('C', [10], 25),('D', [11], 25)]
rows = [('-IPTG (1)', 0, 'r'),\
        ('+IPTG (1)', 2, 'g'),\
        ('-IPTG (2)', 1, 'r--'),\
        ('+IPTG (2)', 3, 'g--')]

for (plot_index, plot_label, y_min, y_max) in plot_types:
    sys.stderr.write("Plotting " + plot_label + " ... \n")
    for i in range(len(columns)):
        fig = figure()
        (name, c_list, max_T) = columns[i]
        
        gr_mat = zeros((len(c_list), len(rows)))
        for j in range(len(c_list)):
            for k in range(len(rows)):
                if (plot_index == "divide"):
                    (time0, values0) = vp.get_data(0, rows[k][1], c_list[j])
                    (time1, values1) = vp.get_data(1, rows[k][1], c_list[j])
                    time = time0
                    values = values1/values0                    
                else:
                    (time, values) = vp.get_data(plot_index, rows[k][1], c_list[j])
                plot(time, values, rows[k][2])
                gr_mat[j, k] = vp.fit_growth(time, values, fit_window_size, fit_start_threshold)
        
        legend_text = []
        for k in range(len(rows)):
            (leg, r, color) = rows[k]
            max_dt = log(2)*60/min(gr_mat[:,r])
            min_dt = log(2)*60/max(gr_mat[:,r])
            legend_text.append("%s (%.0f-%0.f min)" % (leg, min_dt, max_dt))
        legend(legend_text, loc='lower right')
        title(name)
        xlabel('Time (hr)')
        ylabel(plot_label)
        yscale('log')
        axis([0, max_T, y_min, y_max])
        pp.savefig(fig)

pp.close()