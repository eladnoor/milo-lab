import pylab
import random
import matplotlib.pyplot as plt
from scipy import stats


def cdf(v, label, style='b', show_median=False, 
        std=None, y_offset=0, figure=None):
    if figure is None:
        figure = pylab.figure()
        
    new_v_nonan = sorted([x for x in v if x is not None])
    l = len(new_v_nonan)
    if l < 2:
        pylab.plot([0, 0], [0, 0], label=label, figure=figure)
    else:
        y = [float(x)/(float(l)-1.0) for x in range(l)]
        pylab.plot(new_v_nonan, y, style, label=label, figure=figure)
        if std is not None:
            m = pylab.median(new_v_nonan)
            pylab.plot([m-std, m+std], [0.4999 + y_offset, 0.5001 + y_offset], style, lw=1.5, ls='-', figure=figure)
            pylab.plot([m-std-0.0001, m-std+0.0001], [0.485 + y_offset, 0.515 + y_offset], style, lw=1.5, ls='-', figure=figure)
            pylab.plot([m+std-0.0001, m+std+0.0001], [0.485 + y_offset, 0.515 + y_offset], style, lw=1.5, ls='-', figure=figure)
        if show_median:
            m = pylab.median(new_v_nonan)
            pylab.plot([m, m], [0, 0.5], style, figure=figure)
            
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

def binned_plot(x, y, bins, y_type='mean', figure=None, plot_counts=True):
    bins_array = pylab.array([min(x)-1e-14] + list(sorted(bins)) + [max(x)-1e-14])
    binned_y = {}
    for i in xrange(len(x)):
        bin_index = pylab.find(bins_array < x[i]).max()
        binned_y.setdefault(bin_index, []).append(y[i])
    
    y_count = []
    y_vec = []
    for j in xrange(len(bins) + 1):
        if j in binned_y:
            binned_y[j] = pylab.array(binned_y[j])
            y_count.append(len(binned_y[j]))
            if y_type == 'mean':
                y_vec.append(pylab.mean(binned_y[j]))
            elif y_type == 'rmse':
                y_vec.append(pylab.rms_flat(binned_y[j]))
            elif y_type == 'std':
                y_vec.append(pylab.std(binned_y[j]))
        else:
            y_count.append(0)
            y_vec.append(0.0)
    
    bin_width = bins_array[1:] - bins_array[0:-1]
    bin_center = (bins_array[1:] + bins_array[0:-1])/2
    
    if not figure:
        figure = pylab.figure()
    figure.hold(True)
    pylab.bar(left=bins_array[0:-1], height=y_vec, width=bin_width, figure=figure)
    for i in xrange(len(bins) + 1):
        if y_count[i] > 0:
            pylab.text(bin_center[i], y_vec[i], '%d' % y_count[i], horizontalalignment='center', fontsize='small')
        
if __name__ == "__main__":
    x = [0, 1, 0.5, 1.1, 2.5, 6, 4,  7, 5.5]
    y = [0.1, 0.2, 0.1, 0.6,   1.2,   1.3, 1.1, 1.4, 0.1, 0.5]
    bins = [-0.5, 1.5, 3.5, 5.5, 7.5]
    fig = pylab.figure()
    binned_plot(x, y, bins, y_type='rmse', figure=fig)
    pylab.show()