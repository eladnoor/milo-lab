#!/usr/bin/python

COLORS = 'grbcmykw'

def GetColor(i):
	index = i % len(COLORS)
	return COLORS[index]


def ColorMap(items):
	return dict((item, GetColor(i)) for i, item in enumerate(items))
