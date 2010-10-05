from victor_parser import VictorParser
from pylab import *
from matplotlib import font_manager
from matplotlib.backends.backend_pdf import PdfPages
import sys, os

if (not os.path.exists('../res')):
    os.mkdir('../res')

vp = VictorParser()
fname = "../data/OD_eGFP_20100711_190.xls"
vp.parse_excel(fname)

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

columns = [('4K5',      0, 'y'),\
           ('GFP-A',     1, 'r'),\
           ('GFP-B',       2, 'g'),\
           ('GFP-C', 3, 'c'),\
           ('GFP-D',       4, 'c:'),\
           ('BLANK',  5, 'b')]

pp = PdfPages('../res/2010-07-11_lior.pdf')

rows = [('-IPTG', [0,1], 30), ('+IPTG', [3,4], 30)]
for i in range(len(rows)):
    fig = figure()
    (name, r_list, max_T) = rows[i]
    
    for j in range(len(r_list)):
        for k in range(len(columns)):
            (t, od) = vp.get_data(0, r_list[j], columns[k][1])
            (t, gfp) = vp.get_data(1, r_list[j], columns[k][1])
            plot(t, gfp/od, columns[k][2])
    
    legend_text = []
    for k in range(len(columns)):
        (leg, r, color) = columns[k]
        legend_text.append("%s" % (leg))
    legend(legend_text, loc='upper right')
    title(name)
    xlabel('Time (hr)')
    ylabel('GFP/OD')
    yscale('log')
    #axis([0, max_T, 5e3, 1e5])
    pp.savefig(fig)

pp.close()