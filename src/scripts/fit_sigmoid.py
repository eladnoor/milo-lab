import numpy as np
import pylab
import scipy.optimize


def sigmoid(p,x):
    x0,y0,c,k=p
    y = c / (1 + np.exp(-k*(x-x0))) + y0
    return y


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

def growth_rate(params):
    x0,y0,c,k = params
    ymid = y0 + (c / 2.0)
    xmid = (-1/k)*pylab.log((c/(ymid-y0)) - 1) + x0 
    e = pylab.exp(-k*(xmid-x0))
    rate = -k*e / (1+e)**2
    return xmid, ymid, rate
    
def stationary_level(params):
    x0,y0,c,k = params
    return c + y0

def yoffset(params):
    x0,y0,c,k = params
    return y0

p = (0.0, 1.0, 10, 1.0)
high_rate = pylab.arange(-10, 10, 0.25)
low_rate = pylab.arange(-10, 10, 0.5)
vals = sigmoid(p, high_rate)
noisy_vals = noisy_sigmoid(p, low_rate)

pylab.subplot(211)
pylab.plot(high_rate, vals, label='Sigmoid')
pylab.plot(low_rate, noisy_vals, 'r*', label='Samples')

p = fit_sigmoid(low_rate, noisy_vals)
stationary = pylab.zeros(len(high_rate)) + stationary_level(p)
pylab.plot(high_rate, stationary, 'k--', label="Fit stationary")

yoff = pylab.zeros(len(high_rate)) + yoffset(p)
pylab.plot(high_rate, yoff, 'm--', label="Fit Y offsett")

xmid, ymid, slope = growth_rate(p)
pylab.plot([xmid], [ymid], 'mo')
pxp=sigmoid(p, high_rate)
pylab.plot(high_rate, pxp, 'g', label='Fit')

pylab.subplot(212)
residuals = pylab.absolute(noisy_vals - sigmoid(p, low_rate))
average_resid = pylab.mean(residuals)
print average_resid
pylab.stem(range(len(residuals)), residuals)
pylab.show()
