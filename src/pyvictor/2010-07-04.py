from victor_parser import VictorParser
from pylab import *
from matplotlib import font_manager
from matplotlib.backends.backend_pdf import PdfPages
import sys

vp = VictorParser()
fname = "../data/Elad's OD600_20100701_188.xls"
vp.parse_excel(fname)
#(t, v) = vp.get_data(0, 0, 0)
#fit_growth(t, v, 1.5, plot_figure=True)
#sys.exit(0)

#rcParams['text.usetex'] = True
rcParams['legend.fontsize'] = 8
#rcParams['font.family'] = 'sans-serif'
#rcParams['font.size'] = 8
#rcParams['lines.linewidth'] = 0.3
#rcParams['lines.markersize'] = 2
#rcParams['figure.figsize'] = [5, 10]
#rcParams['figure.subplot.hspace'] = 0.3
#figure()

fit_window_size = 1.5 # hours
fit_start_threshold = 0.01

rows = [('SOC',      0, 'y'),\
        ('casa',     1, 'r'),\
        ('DO',       2, 'g'),\
        ('AA+dntos', 3, 'c'),\
        ('AA',       4, 'c:'),\
        ('acetate',  5, 'b'),\
        ('glucose',  6, 'm'),\
        ('glucose stationary', 7, 'm:')]

pp = PdfPages('../res/2010-07-04.pdf')

columns = [('PRK', [0,1,2], 15), ('GFP', [3,4,5], 20), ('control', [6,7], 50)]
for i in range(len(columns)):
    fig = figure()
    (name, c_list, max_T) = columns[i]
    
    gr_mat = zeros((len(c_list), len(rows)))
    for j in range(len(c_list)):
        for k in range(len(rows)):
            (t, v) = vp.get_data(0, rows[k][1], c_list[j])
            plot(t, v-min(v), rows[k][2])
            gr_mat[j, k] = vp.fit_growth(t, v, fit_window_size, fit_start_threshold)
            #gr_mat[j, k] = vp.fit_growth2(t, v)
    
    legend_text = []
    for k in range(len(rows)):
        (leg, r, color) = rows[k]
        max_dt = log(2)*60/min(gr_mat[:,r])
        min_dt = log(2)*60/max(gr_mat[:,r])
        legend_text.append("%s (%.0f-%0.f min)" % (leg, min_dt, max_dt))
    legend(legend_text, loc='lower right')
    title(name)
    xlabel('Time (hr)')
    ylabel('GFP')
    yscale('log')
    axis([0, max_T, 0.001, 1])
    pp.savefig(fig)

pp.close()