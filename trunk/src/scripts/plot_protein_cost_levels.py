import numpy as np
import pylab

from matplotlib.ticker import MultipleLocator, LogFormatterMathtext, FormatStrFormatter

vtotal = 5.0 / (60**2 * 1000)
factor = vtotal * 100.0 * 100.0 / 55.0
min_ed_vals = np.array([2212, 3052]) * factor
min_emp_vals = np.array([2217, 5854]) * factor
ed_vals = np.array([3516, 4978]) * factor
emp_vals = np.array([1.197e4, 2.589e4]) * factor
left = np.array(range(2))
xticks = ['Avg', 'Kin']

ytick_formatter = FormatStrFormatter('%d%%')
f = pylab.figure()

width=0.35
offset=0.15
pylab.bar(left+offset, ed_vals, color='m', width=width, edgecolor='w')
pylab.bar(left+offset, min_ed_vals, color='#cccccc', width=width, edgecolor='w')
pylab.bar(left+offset+width, emp_vals, color='g', width=width, edgecolor='w')
pylab.bar(left+offset+width, min_emp_vals, color='#cccccc', width=width, edgecolor='w')
#pylab.ylim((100, 3e4))
pylab.xticks(left + 0.5, xticks)
ax = f.get_axes()
ax[0].yaxis.set_major_formatter(ytick_formatter)

pylab.show()