#!/usr/bin/python

COLORS = 'grbcmyk'

import colorsys
import numpy as np
import matplotlib.pyplot as plt

def RgbToHex(rgb):
	r, g, b = rgb
	r *= 255
	g *= 255
	b *= 255
	return '#%02x%02x%02x' % (r,g,b)
	

def ColorMap(items, saturation=0.7, value=0.95, hues=None):
	if hues is None:
		n = len(items)
		hues = np.arange(float(n)) / float(n)
	f = lambda h: RgbToHex(colorsys.hsv_to_rgb(h, saturation, value)) 
	rgbs = map(f, hues)
	return dict(zip(items, rgbs))

def LayeredColorMap(items, n_minor, saturation=(0.4, 0.7), value=(0.95, 0.5), 
				    hues=None):
	n = len(saturation)
	colormaps = [ColorMap(items, saturation[i], value[i], hues) for i in xrange(n)]
	colormap = {}
	for item in items:
		colormap[item] = [colormaps[i % n][item] for i in xrange(n_minor)]
	return colormap

def test():
	items = ['A', 'B', 'C', 'D', 'E']
	hues = [0.05, 0.22, 0.5, 0.7, 0.9]
	n = 5
	m = 10
	colormap = LayeredColorMap(items, n, saturation=(0.4, 1.0), value=(1.0, 0.4), hues=hues)
	data = abs(np.random.randn(len(items)*n, m))
	bottom = np.zeros((1, m))

	ind = np.arange(m)    # the x locations for the groups
	width = 0.35          # the width of the bars
	
	for i, item in enumerate(items):
		for j in xrange(n):
			plt.bar(ind, data[n*i+j, :], width, color=colormap[item][j],
				    label=item + '%d' % j, bottom=bottom[0, :], edgecolor='none')
			bottom += data[n*i+j, :]

	top = max(bottom[0, :]) - bottom
	plt.bar(ind, top[0, :], width, color='#AAAAAA',
		    label=item + '%d' % j, bottom=bottom[0, :], edgecolor='none')
	
	plt.legend()
	plt.show()

if __name__ == "__main__":
	test()