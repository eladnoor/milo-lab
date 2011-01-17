import pylab

def cdf(v, label, style='b', show_median=False):
    new_v = []
    for x in v:
        if (x != None):
            new_v.append(x)
    l = len(new_v)
    if (l < 2):
        pylab.plot([0, 0], [0, 0], label=label)
    else:
        new_v.sort()
        
        m = pylab.median(v)
        pylab.plot(new_v, [float(x)/(float(l)-1) for x in range(l)], style, label=label)
        if show_median:
            pylab.hold(True)
            pylab.plot([m, m], [0, 0.5], style)
    
