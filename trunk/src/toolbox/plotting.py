import random
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.patches as patch
from matplotlib.mlab import rms_flat

def cdf(v, label, style='b', show_median=False, 
        std=None, y_offset=0, figure=None):
    if figure is None:
        figure = plt.figure()
        
    new_v_nonan = sorted([x for x in v if x is not None])
    l = len(new_v_nonan)
    if l < 2:
        plt.plot([0, 0], [0, 0], label=label, figure=figure)
    else:
        y = [float(x)/(float(l)-1.0) for x in range(l)]
        plt.plot(new_v_nonan, y, style, label=label, figure=figure)
        if std is not None:
            m = np.median(new_v_nonan)
            plt.plot([m-std, m+std], [0.4999 + y_offset, 0.5001 + y_offset], style, lw=1.5, ls='-', figure=figure)
            plt.plot([m-std-0.0001, m-std+0.0001], [0.485 + y_offset, 0.515 + y_offset], style, lw=1.5, ls='-', figure=figure)
            plt.plot([m+std-0.0001, m+std+0.0001], [0.485 + y_offset, 0.515 + y_offset], style, lw=1.5, ls='-', figure=figure)
        if show_median:
            m = np.median(new_v_nonan)
            plt.plot([m, m], [0, 0.5], style, figure=figure)
            
# Yaniv: Not used  
def bootstrap(values_map,colors_map, func=np.median, reps=10000, sample_frac=1, name='Median'):

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
        sample_mean = np.mean(sample_vals)
        sample_std = np.std(sample_vals)
        sample_label = r'$%s:\ \mu=%.2f,\ \sigma=%.2f$' % (str(label), sample_mean, sample_std)
            
        plt.hist(sample_vals, 100, normed=1, ec=None, facecolor=colors_map.get(label, 'b'), alpha=0.75, 
                 label=sample_label)
        plt.plot([m, m], [0, 0.5], 'b--') #colors_map.get(label, 'b'))
    plt.xlabel(name)
    plt.ylabel('Fraction')
    
    plt.legend(loc='upper left')

def binned_plot(x, y, bins, y_type='mean', figure=None, plot_counts=True):
    bins_array = np.array([min(x)-1e-14] + list(sorted(bins)) + [max(x)-1e-14])
    binned_y = {}
    for i in xrange(len(x)):
        bin_index = max(np.nonzero(bins_array < x[i])[0])
        binned_y.setdefault(bin_index, []).append(y[i])
    
    y_count = []
    y_vec = []
    for j in xrange(len(bins) + 1):
        if j in binned_y:
            binned_y[j] = np.array(binned_y[j])
            y_count.append(len(binned_y[j]))
            if y_type == 'mean':
                y_vec.append(np.mean(binned_y[j]))
            elif y_type == 'rmse':
                y_vec.append(rms_flat(binned_y[j]))
            elif y_type == 'std':
                y_vec.append(np.std(binned_y[j]))
        else:
            y_count.append(0)
            y_vec.append(0.0)
    
    bin_width = bins_array[1:] - bins_array[0:-1]
    bin_center = (bins_array[1:] + bins_array[0:-1])/2
    
    if not figure:
        figure = plt.figure()
    figure.hold(True)
    plt.bar(left=bins_array[0:-1], height=y_vec, width=bin_width, figure=figure)
    for i in xrange(len(bins) + 1):
        if y_count[i] > 0:
            plt.text(bin_center[i], y_vec[i], '%d' % y_count[i], horizontalalignment='center', fontsize='small')
        

def bubble_plot(x, y, s, e=None, c=None, figure=None):
    """
        Inputs:
            x - a list of x-values for the bubbles
            y - a list of y-values for the bubbles
            s - a list of the areas of the bubbles
            e - a list of pairs of min and max areas for the error bars
            c - a list of RGB-tuples
            
    """
    if figure is None:
        figure = plt.figure()
    figure.clf()
    ax = fig.add_subplot(111, autoscale_on=True)
    figure.hold(True)
    N = len(x)
    
    x_edge = np.zeros((N))
    if e:
        x_err = np.zeros((2, N))
    for i in xrange(N):
        radius = np.sqrt(s[i]/np.pi)
        x_edge[i] = x[i] + radius
        if c:
            color = c[i]
        else:
            color = (0, 0, 1)
        circ = patch.Circle((x[i], y[i]), radius=radius, figure=figure,
                            linewidth=0, facecolor=color)
        ax.add_patch(circ)
        if e:
            x_err[0, i] = radius - np.sqrt(e[i][0]/np.pi)
            x_err[1, i] = np.sqrt(e[i][1]/np.pi) - radius
    
    if e:
        plt.errorbar(x=x_edge, y=y, xerr=x_err, fmt=None, ecolor='black')
    ax.axes.relim()
    ax.axes.autoscale_view()

if __name__ == "__main__":
    fig = plt.figure(figsize=(6,6))

    #x = [0, 1, 0.5, 1.1, 2.5, 6, 4,  7, 5.5]
    #y = [0.1, 0.2, 0.1, 0.6,   1.2,   1.3, 1.1, 1.4, 0.1, 0.5]
    #bins = [-0.5, 1.5, 3.5, 5.5, 7.5]
    #binned_plot(x, y, bins, y_type='rmse', figure=fig)
    bubble_plot([0, 0, 5, 5], [0, 5, 0, 5], [1, 2, 2, 1], 
                e=[(0.8, 1.2), (1.6, 2.5), (1.8, 2.1), (0.5, 1.5)],
                c=[(0.4, 0.6, 0.0), (0.4, 0.2, 0.0), (0.0, 0.6, 0.5), (0.9, 0.9, 0.0)],
                figure=fig)
    plt.show()