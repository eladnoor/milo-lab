from victor_parser import VictorParser
from pylab import *
from matplotlib import font_manager
from matplotlib.backends.backend_pdf import PdfPages
import sys, os

if (not os.path.exists('../res')):
    os.mkdir('../res')

vp = VictorParser()
name = "RBS-GFP 28-07-10"
fname = "../data/victor/%s.xls" % name
if (not os.path.exists(fname)):
    raise Exception("Cannot locate the Excel file: " + fname)
vp.parse_excel(fname)
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

plot_types = [(0, "OD", 0.001, 1), (1, "GFP", 1e3, 1e6), ("divide", "GFP/OD", 1e4, 2e6)]

plots = [] # (title, victor_index, (t_min, t_max), (y_min, y_max), y_label, 

plots.append(('Blank', 0, (0, 25), (1e-3, 1), 'OD', {'-IPTG':('r', [(0,0),(1,0),(2,0)]), '+IPTG':('g', [(3,0),(4,0),(5,0)])}))
plots.append(('RBS-A', 0, (0, 25), (1e-3, 1), 'OD', {'-IPTG':('r', [(0,1),(1,1),(2,1)]), '+IPTG':('g', [(3,1),(4,1),(5,1)])}))
plots.append(('RBS-B', 0, (0, 25), (1e-3, 1), 'OD', {'-IPTG':('r', [(0,2),(1,2),(2,2)]), '+IPTG':('g', [(3,2),(4,2),(5,2)])}))
plots.append(('RBS-C', 0, (0, 25), (1e-3, 1), 'OD', {'-IPTG':('r', [(0,3),(1,3),(2,3)]), '+IPTG':('g', [(3,3),(4,3),(5,3)])}))
plots.append(('RBS-D', 0, (0, 25), (1e-3, 1), 'OD', {'-IPTG':('r', [(0,4),(1,4),(2,4)]), '+IPTG':('g', [(3,4),(4,4),(5,4)])}))
plots.append(('Plasmid', 0, (0, 25), (1e-3, 1), 'OD', {'-IPTG':('r', [(0,5),(1,5),(2,5)]), '+IPTG':('g', [(3,5),(4,5),(5,5)])}))

plots.append(('Blank', 1, (0, 25), (1e3, 1e6), 'GFP', {'-IPTG':('r', [(0,0),(1,0),(2,0)]), '+IPTG':('g', [(3,0),(4,0),(5,0)])}))
plots.append(('RBS-A', 1, (0, 25), (1e3, 1e6), 'GFP', {'-IPTG':('r', [(0,1),(1,1),(2,1)]), '+IPTG':('g', [(3,1),(4,1),(5,1)])}))
plots.append(('RBS-B', 1, (0, 25), (1e3, 1e6), 'GFP', {'-IPTG':('r', [(0,2),(1,2),(2,2)]), '+IPTG':('g', [(3,2),(4,2),(5,2)])}))
plots.append(('RBS-C', 1, (0, 25), (1e3, 1e6), 'GFP', {'-IPTG':('r', [(0,3),(1,3),(2,3)]), '+IPTG':('g', [(3,3),(4,3),(5,3)])}))
plots.append(('RBS-D', 1, (0, 25), (1e3, 1e6), 'GFP', {'-IPTG':('r', [(0,4),(1,4),(2,4)]), '+IPTG':('g', [(3,4),(4,4),(5,4)])}))
plots.append(('Plasmid', 1, (0, 25), (1e3, 1e6), 'GFP', {'-IPTG':('r', [(0,5),(1,5),(2,5)]), '+IPTG':('g', [(3,5),(4,5),(5,5)])}))

plots.append(('Blank', 'divide', (0, 25), (1e4, 1e7), 'GFP/OD', {'-IPTG':('r', [(0,0),(1,0),(2,0)]), '+IPTG':('g', [(3,0),(4,0),(5,0)])}))
plots.append(('RBS-A', 'divide', (0, 25), (1e4, 1e7), 'GFP/OD', {'-IPTG':('r', [(0,1),(1,1),(2,1)]), '+IPTG':('g', [(3,1),(4,1),(5,1)])}))
plots.append(('RBS-B', 'divide', (0, 25), (1e4, 1e7), 'GFP/OD', {'-IPTG':('r', [(0,2),(1,2),(2,2)]), '+IPTG':('g', [(3,2),(4,2),(5,2)])}))
plots.append(('RBS-C', 'divide', (0, 25), (1e4, 1e7), 'GFP/OD', {'-IPTG':('r', [(0,3),(1,3),(2,3)]), '+IPTG':('g', [(3,3),(4,3),(5,3)])}))
plots.append(('RBS-D', 'divide', (0, 25), (1e4, 1e7), 'GFP/OD', {'-IPTG':('r', [(0,4),(1,4),(2,4)]), '+IPTG':('g', [(3,4),(4,4),(5,4)])}))
plots.append(('Plasmid', 'divide', (0, 25), (1e4, 1e7), 'GFP/OD', {'-IPTG':('r', [(0,5),(1,5),(2,5)]), '+IPTG':('g', [(3,5),(4,5),(5,5)])}))

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

    legend(legend_text, loc='lower right')
    yscale('log')
    axis([t_range[0], t_range[1], y_range[0], y_range[1]])
    pp.savefig(fig)

pp.close()