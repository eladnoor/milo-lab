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


def MakeBinary(a):
	return a.lower() == 'true'

def MergeBinary(a, b, a_val, b_val):
	abin = MakeBinary(a)
	bbin = MakeBinary(b)

	if abin and bbin:
		return 'BOTH'
	if abin:
		return a_val
	if bbin:
		return b_val
	return 'NEITHER'


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
		if filter_values and value in filter_values:
			continue

		all_vals.add(value)
		pairmap.setdefault(a, []).append(value)
	
	all_weights = {}
	all_counts = {}
	for k, v in pairmap.iteritems():
		counter = Counter(v)
		all_counts[k] = counter
		
		total = float(sum(counter.values()))
		weights = {}
		for k2, v2 in counter.iteritems():
			weights[k2] = float(v2) / total
		all_weights[k] = weights
		
		
	weight_array = []
	count_array = []
	for a_val in all_a:
		weight_array.append(all_weights.get(a_val, {}))
		count_array.append(all_counts.get(a_val, {}))

	# Plot circle scatter.
	axes = pylab.axes()

	colormap = ColorMap(all_vals)
	ind = pylab.arange(len(all_a))
	current_bottom = pylab.zeros(len(all_a))
	for key in sorted(all_vals):
		heights = pylab.array([w.get(key, 0.0) for w in weight_array])
		
		pylab.bar(ind, heights, color=colormap[key],
				  bottom=current_bottom, label=key, width=0.5)
		
		counts = [c.get(key, 0) for c in count_array]
		for x, y, count, weight in zip(ind, current_bottom, counts, heights):
			if count == 0:
				continue
			
			txt = '%d; %.1f%%' % (count, 100*weight)
			pylab.text(x + 0.25, y + 0.02, txt, ha='center')
			
		current_bottom += heights

	title = options.title or '%s vs. %s' % (options.first_col, options.second_col)
	pylab.title(title)

	size_12 = FontProperties(size=12)
	a_labels = [a or "None given" for a in all_a]
	pylab.xticks(ind + 0.25, a_labels, fontproperties=size_12)
	axes.yaxis.set_major_locator(pylab.NullLocator())

	pylab.legend()
	pylab.show()	


if __name__ == '__main__':
	Main()
