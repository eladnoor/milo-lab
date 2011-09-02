import numpy as np
import pylab
import scipy.optimize
import sys


def sigmoid(p,x):
    x0,y0,c,k=p
    y = c / (1 + np.exp(-k*(x-x0))) + y0
    return y

def linear(m, x, b):
    return m*x+b

def noisy_sigmoid(p, x):
    r = sigmoid(p,x)
    noise = 0.2*pylab.randn(r.size)
    return r + noise


def residuals(p,x,y):
    return y - sigmoid(p,x)


def fit_sigmoid(x, y):
    p_guess=(np.median(x), np.median(y),
             max(y), min(y))
    p, cov, infodict, mesg, ier = scipy.optimize.leastsq(
            residuals, p_guess, args=(x,y),
            full_output=1)
    return p


def midpoint(params):
    x0,y0,c,k = params
    ymid = y0 + (c / 2.0)
    xmid = (-1/k)*pylab.log((c/(ymid-y0)) - 1) + x0
    return xmid, ymid 

def growth_rate(params):
    x0,y0,c,k = params
    ymid, ymid = midpoint(params) 
    e = pylab.exp(-k*(xmid-x0))
    rate = k*c*e / (1+e)**2
    offset = ymid - rate*xmid
    return rate, offset
    
def stationary_level(params):
    x0,y0,c,k = params
    return c + y0

def yoffset(params):
    x0,y0,c,k = params
    return y0


noisy_vals = pylab.array(
[4.62e-2,
5.54e-2,
5.7e-2, 5.49e-2, 5.93e-2, 6.21e-2, 6.82e-2, 7.75e-2, 9.13e-2, 0.1091, 0.127, 0.1481, 0.17, 0.1989, 0.2339, 0.2719, 0.3174, 0.3489, 0.383, 0.417, 0.4488, 0.4889, 0.5332, 0.5679, 0.6023, 0.6305, 0.6551, 0.6795, 0.7075, 0.7318, 0.7591, 0.7798, 0.8016, 0.8254, 0.8434, 0.8646, 0.8816, 0.9015, 0.9159, 0.9296, 0.9356, 0.9373, 0.9364, 0.9375, 0.9368, 0.9368, 0.9385, 0.9395, 0.941, 0.9376, 0.9355, 0.9361, 0.9345, 0.937, 0.9369, 0.9329, 0.9331, 0.9335, 0.9336, 0.9331, 0.9308, 0.9318, 0.93, 0.9296, 0.9315, 0.9283, 0.9285, 0.9273, 0.9275, 0.9274, 0.9275, 0.9248])

noisy_vals = pylab.log(noisy_vals)
"""
noisy_vals = pylab.array([0.0314,
0.032,
0.0397,
0.0408,
0.0412,
0.0415,
0.0421,
0.0426,
0.0434,
0.0443,
0.0456,
0.0472,
0.0477,
0.0506,
0.0549,
0.0579,
0.064,
0.0704,
0.079,
0.0921,
0.1101,
0.1334,
0.1675,
0.2036,
0.2493,
0.2942,
0.3393])
"""
times = pylab.arange(len(noisy_vals))


#p = (0.0, 1.0, 10, 1.0)
high_rate = pylab.arange(min(times), 1.5*max(times), 0.25)
low_rate = pylab.arange(min(times), 1.5*max(times), 0.5)
#vals = sigmoid(p, high_rate)
#noisy_vals = noisy_sigmoid(p, low_rate)

pylab.subplot(211)
#pylab.plot(high_rate, vals, label='Sigmoid')
pylab.plot(times, noisy_vals, 'r*', label='Samples')

p = fit_sigmoid(times, noisy_vals)
stationary = pylab.zeros(len(times)) + stationary_level(p)
pylab.plot(times, stationary, 'k--', label="Fit stationary")

yoff = pylab.zeros(len(times)) + yoffset(p)
print yoff
pylab.plot(times, yoff, 'm--', label="Fit Y offsett")

xmid, ymid = midpoint(p)
slope, offset = growth_rate(p)
print slope, offset
pylab.plot([xmid], [ymid], 'mo', label='Inflection point')


l = linear(slope, high_rate, offset)
growth_x = []
growth_y = []
for x, y in zip(high_rate, l):
    if y < yoff[0] or y > stationary[0]:
        continue
    growth_x.append(x)
    growth_y.append(y)
        
pylab.plot(growth_x, growth_y, 'k-', label='Growth phase')
pxp=sigmoid(p, high_rate)
pylab.plot(high_rate, pxp, 'g', label='Fit')
pylab.legend(loc="lower right")

pylab.subplot(212)
residuals = pylab.absolute(noisy_vals - sigmoid(p, times))
average_resid = pylab.mean(residuals)
pylab.stem(range(len(residuals)), residuals)
pylab.show()
