from victor_parser import VictorParser
from pylab import *
from matplotlib import font_manager
import sys

vp = VictorParser()
fname = "../data/Elad's OD600_20100629_186.xls"
vp.parse_excel(fname)
#(t, v) = vp.get_data(0, 0, 0)
#fit_growth(t, v, 1.5, plot_figure=True)
#sys.exit(0)

figure(figsize=(5, 7))
rcParams['text.usetex'] = True
rcParams['legend.fontsize'] = 6
rcParams['font.family'] = 'sans-serif'
rcParams['font.size'] = 8
rcParams['lines.linewidth'] = 0.6
rcParams['lines.markersize'] = 2

fit_window_size = 1.5 # hours
fit_start_threshold = 0.01
leg = ['SOC+kana', 'M9+CAS+Glu+kana', 'M9+AA+Nt+Glu+kana', 'M9+Glu+kana', 'M9+Acetate+kana']
color = ['g','r','b','m','c']
gr_mat = zeros((5,10))

subplot(2,1,1)
for c in range(0, 5):
    for r in range(5):
        (t, v) = vp.get_data(0, r, c)
        plot(t, v-min(v), color[r])
        (gr, err) = vp.fit_growth(t, v, fit_window_size, fit_start_threshold)
        gr_mat[r,c] = gr
legend([leg[r] + " (%.0f min)" % (log(2)*60/mean(gr_mat[r,0:5])) for r in range(5)], loc='lower right')
title('PRK')
xlabel('Time (hr)')
ylabel('OD')
yscale('log')
axis([0, t.max(), 0.001, 1])


subplot(2,1,2)
for c in range(5, 10):
    for r in range(5):
        (t, v) = vp.get_data(0, r, c)
        plot(t, v-min(v), color[r])
        (gr, err) = vp.fit_growth(t, v, fit_window_size, fit_start_threshold)
        gr_mat[r,c] = gr
legend([leg[r] + " (%.0f min)" % (log(2)*60/mean(gr_mat[r,5:10])) for r in range(5)], loc='lower right')
yscale('log')
title('GFP')
xlabel('Time (hr)')
ylabel('OD')
axis([0, t.max(), 0.001, 1])
savefig('../res/2010-06-29.pdf', format='pdf')