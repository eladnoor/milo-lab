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


noisy_vals =[0.0374,0.0363,0.0396,0.0396,0.0394,0.0393,0.0394,0.0394,0.0395,0.0395,
             0.0397,0.0397,0.0398,0.0395,0.0399,0.0399,0.04,0.0402,0.04,0.0401,0.0402,
             0.0405,0.0403,0.0406,0.0402,0.0407,0.0409,0.041,0.0411,0.0412,0.0414,0.0418,
             0.0417,0.0421,0.0422,0.043,0.0429,0.0438,0.044,0.0445,0.0446,0.0453,0.0461,0.0477,
             0.0485,0.0494,0.0504,0.0519,0.0537,0.0564,0.0577,0.061,0.0645,0.0684,
             0.0732,0.0791,0.0851,0.0915,0.1012,0.1102,0.1212,0.1333,0.1455,0.1534,
             0.1615,0.1762,0.1848,0.1962,0.2133,0.2276,0.2386,0.2501,0.2645,0.273,
             0.2845,0.2947,0.3062,0.3115,0.3189,0.3274,0.3355,0.3401,0.3458,0.3464,
             0.3494,0.358,0.3508,0.3537,0.3406,0.334,0.3411,0.3301,0.326]

noisy_vals = pylab.log(noisy_vals)
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
