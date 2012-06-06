import numpy as np
import pylab

from matplotlib.ticker import MultipleLocator, LogFormatterMathtext, FormatStrFormatter


ed_levels_avg  = [0.082900816606927483, 0.089511793895029848, 0.083309766755565826, 0.084484853574007535, 0.14724712001615858, 0.14521870622850469, 0.16355192943011002, 0.12294765262467099, 0.13325646618501297, 0.09301821831723682, 0.27515909697097701]
emp_levels_avg = [0.082900816606927483, 0.0871928878523133, 0.085654802175554628, 0.34537367640444849, 0.35277365823506412, 0.68120320079172325, 1.6172738321832865, 0.48799836557984327, 0.71534993991045293, 0.2485110415876649, 0.17902230256282661]

ed_sum  = sum(ed_levels_avg)
emp_sum = sum(emp_levels_avg)
print 'EMP/ED', (emp_sum / ed_sum)

def format_tick(val, pos):
    if val <= 0:
        return val
    exp = int(np.floor(np.log10(val)))
    const = val / 10**exp
    return '$%.1f\\times10^%d$' % (const, exp)

ytick_formatter = FormatStrFormatter('%d%%')
ytick_locator   = MultipleLocator(100)


f = pylab.figure()
ax = f.get_axes()
pylab.plot(np.cumsum(ed_levels_avg), 'm-', lw=3, label='ED')
pylab.plot(np.cumsum(emp_levels_avg), 'g-', lw=3, label='EMP')
pylab.legend()
pylab.xticks(pylab.arange(0, len(emp_levels_avg)),
             pylab.arange(1, len(emp_levels_avg) + 1))
#ax[0].yaxis.set_major_locator(ytick_locator)
ax[0].yaxis.set_major_formatter(ytick_formatter)

pylab.show()

