import numpy as np
import pylab
import scipy.optimize


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
    p_guess=(np.median(low_rate), np.median(noisy_vals),
             max(noisy_vals), min(noisy_vals))
    p, cov, infodict, mesg, ier = scipy.optimize.leastsq(
            residuals, p_guess, args=(low_rate,noisy_vals),
            full_output=1, warning=True)
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
times = pylab.arange(len(noisy_vals))


#p = (0.0, 1.0, 10, 1.0)
#high_rate = pylab.arange(-10, 10, 0.25)
#low_rate = pylab.arange(-10, 10, 0.5)
#vals = sigmoid(p, high_rate)
#noisy_vals = noisy_sigmoid(p, low_rate)

pylab.subplot(211)
#pylab.plot(high_rate, vals, label='Sigmoid')
pylab.plot(times, noisy_vals, 'r*', label='Samples')

p = fit_sigmoid(times, noisy_vals)
stationary = pylab.zeros(len(times)) + stationary_level(p)
pylab.plot(times, stationary, 'k--', label="Fit stationary")

yoff = pylab.zeros(len(times)) + yoffset(p)
pylab.plot(times, yoff, 'm--', label="Fit Y offsett")

xmid, ymid = midpoint(p)
slope, offset = growth_rate(p)
pylab.plot([xmid], [ymid], 'mo', label='Inflection point')
pylab.plot(high_rate, linear(slope, high_rate, offset), 'k-', label='Growth phase')
pxp=sigmoid(p, high_rate)
pylab.plot(high_rate, pxp, 'g', label='Fit')
pylab.legend(loc="upper left")


pylab.subplot(212)
residuals = pylab.absolute(noisy_vals - sigmoid(p, low_rate))
average_resid = pylab.mean(residuals)
pylab.stem(range(len(residuals)), residuals)
pylab.show()
