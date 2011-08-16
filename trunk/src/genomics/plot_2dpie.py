#!/usr/bin/python

import csv
import itertools
import pylab
import sys

from collections import Counter
from matplotlib.font_manager import FontProperties
from matplotlib.patches import Wedge
from optparse import OptionParser
from colormap import ColorMap


def MakeOpts():
	"""Returns an OptionParser object with all the default options."""
	opt_parser = OptionParser()
	opt_parser.add_option("-i", "--input_filename",
						  dest="input_filename",
						  help="input CSV")
	opt_parser.add_option("-a", "--first_col",
						  dest="first_col",
						  help="First column in 2-way histogram")
	opt_parser.add_option("-b", "--second_col",
						  dest="second_col",
						  help="Second column in 2-way histogram")
	opt_parser.add_option("-c", "--third_col",
						  dest="third_col",
						  help="Third column in 2-way histogram")
	opt_parser.add_option("-v", "--second_col_key",
						  dest="second_col_key",
						  help="First column in 2-way histogram")
	opt_parser.add_option("-w", "--third_col_key",
						  dest="third_col_key",
						  help="First column in 2-way histogram")	
	opt_parser.add_option("-f", "--filter_values",
						  dest="filter_values",
						  help="Value to ignore.")	
	opt_parser.add_option("-t", "--title",
						  dest="title",
						  help="Figure Title")	
	return opt_parser


def PieArray(axes, x_array, y_loc, radius, counts_array, colormap, **kwargs):
    for x, counts in zip(x_array, counts_array):
		total_count = sum(counts.itervalues())
		axes.text(x, y_loc + radius + 0.05, str(total_count), ha='center')

		theta_1 = 0
		for i, (key, value) in enumerate(counts.iteritems()):
			color = colormap[key]
			sweep = (float(value) / float(total_count)) * 360.0
			w = Wedge((x,y_loc), radius, theta_1, theta_1 + sweep, color=color, label=key)
			theta_1 += sweep
			axes.add_patch(w)

    return True


def MakeBinary(a):
	return a.lower() == 'true'

def MergeBinary(a, b, a_val, b_val):
	abin = MakeBinary(a)
	bbin = MakeBinary(b)

	if abin and bbin:
		return [a_val, b_val]
	if abin:
		return [a_val]
	if bbin:
		return [b_val]
	return ['NEITHER']


def Main():
	options, _ = MakeOpts().parse_args(sys.argv)
	assert options.input_filename
	assert options.first_col and options.second_col
	print 'Reading species list from', options.input_filename

	# Read data
	r = csv.DictReader(open(options.input_filename))
	first_col, second_col = options.first_col, options.second_col
	third_col = options.third_col
	second_col_key = options.second_col_key or second_col
	third_col_key = options.third_col_key or third_col
	filter_values = []
	if options.filter_values:
		filter_values = set(map(str.strip, options.filter_values.split(',')))

	pairmap = {}
	all_a = set()
	all_b = set()
	all_vals = set()
	for row in r:
		a, b = row[first_col].strip(), row[second_col].strip()
		c = row[third_col].strip()
		if not a or not b or not c:
			continue
		all_a.add(a)
		all_b.add(b)

		value = MergeBinary(b, c, second_col_key, third_col_key)
		#if filter_values and value in filter_values:
		#	continue

		all_vals.update(value)
		pairmap.setdefault(a, []).extend(value)
	
	a_to_num = dict((v,i) for i,v in enumerate(all_a))

	counts = {}
	for k, v in pairmap.iteritems():
		counter = Counter(v)
		counts[k] = counter

	x_vals = []
	count_array = []
	for a_val in all_a:
		x_vals.append(a_to_num[a_val])
		count_array.append(counts.get(a_val, {}))

	# Plot circle scatter.
	axes = pylab.axes()

	colormap = ColorMap(all_vals)
	y_loc, radius = 1.0, 0.45
	PieArray(axes, x_vals, y_loc=y_loc, radius=radius,
			 counts_array=count_array, colormap=colormap)

	title = options.title or '%s vs. %s' % (options.first_col, options.second_col)
	pylab.title(title)

	size_10 = FontProperties(size=10)
	a_labels = [a or "None given" for a in all_a]
	axes.yaxis.set_major_locator(pylab.NullLocator())
	axes.xaxis.set_major_locator(pylab.NullLocator())
	for x, label in enumerate(a_labels):
		axes.text(x, y_loc - radius - 0.1, str(label), ha='center',
				  va='baseline')

	handles, labels = axes.get_legend_handles_labels()
	mapped_labels = dict(zip(labels, handles))
	labels = sorted(mapped_labels.keys())
	handles = [mapped_labels[k] for k in labels]
	pylab.legend(handles, labels)

	pylab.axis('scaled')
	pylab.axis([-1, len(all_a), 0.0, 2])
	pylab.show()	


if __name__ == '__main__':
	Main()
