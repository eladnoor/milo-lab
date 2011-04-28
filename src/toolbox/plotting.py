import pylab
import random
import matplotlib.pyplot as plt
from scipy import stats


def cdf(v, label, style='b', show_median=False, std=None, y_offset=0):
    new_v = []
    for x in v:
        if (x != None):
            new_v.append(x)
    l = len(new_v)
    if (l < 2):
        pylab.plot([0, 0], [0, 0], label=label)
    else:
        new_v.sort()
                           
        m = stats.cmedian(v)
        y = [float(x)/(float(l)-1) for x in range(l)]
        pylab.plot(new_v, y, style, label=label)
        #pylab.errorbar(new_v, y, xerr=std, style, label=label)
        if std != None:
            pylab.plot([m-std, m+std], [0.4999 + y_offset, 0.5001 + y_offset], style, lw=1.5, ls='-')
            pylab.plot([m-std-0.0001, m-std+0.0001], [0.485 + y_offset, 0.515 + y_offset], style, lw=1.5, ls='-')
            pylab.plot([m+std-0.0001, m+std+0.0001], [0.485 + y_offset, 0.515 + y_offset], style, lw=1.5, ls='-')
        if show_median:
            #pylab.hold(True)
            pylab.plot([m, m], [0, 0.5], style)
            
# Yaniv: Not used  
def bootstrap(values_map,colors_map, func=pylab.median, reps=10000, sample_frac=1, name='Median'):

    for label,vals in values_map.iteritems():
        
        if len(vals) < 20:
            continue
        
        m = func(vals)
        sample_vals = []
        n_to_sample = sample_frac * len(vals)

        i = 0
        while i < reps:
            samples = []
            while len(samples) < n_to_sample:
                samples.append(random.choice(vals))
            sample_vals.append(func(samples))
            i += 1
        sample_mean = pylab.mean(sample_vals)
        sample_std = pylab.std(sample_vals)
        sample_label = r'$%s:\ \mu=%.2f,\ \sigma=%.2f$' % (str(label), sample_mean, sample_std)
            
        plt.hist(sample_vals, 100, normed=1, ec=None, facecolor=colors_map.get(label, 'b'), alpha=0.75, 
                 label=sample_label)
        pylab.plot([m, m], [0, 0.5], 'b--') #colors_map.get(label, 'b'))
    pylab.xlabel(name)
    pylab.ylabel('Fraction')
    
    pylab.legend(loc='upper left')
