#!/usr/bin/python

COLORS = 'grbcmyk'

import colorsys
import numpy

def RgbToHex(rgb):
	r, g, b = rgb
	r *= 255
	g *= 255
	b *= 255
	return '#%02x%02x%02x' % (r,g,b)
	

def ColorMap(items, saturation=0.7, value=0.95):
	n = len(items)
	hues = numpy.arange(float(n)) / float(n)
	f = lambda h: RgbToHex(colorsys.hsv_to_rgb(h, saturation, value)) 
	rgbs = map(f, hues)
	return dict(zip(items, rgbs))
