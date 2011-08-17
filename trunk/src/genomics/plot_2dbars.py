#!/usr/bin/python

import csv
import itertools
import logging
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
	opt_parser.add_option("-o", "--output_filename",
						  dest="output_filename",
						  help="output CSV")
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
	opt_parser.add_option("-x", "--filter_cols",
						  dest="filter_cols",
						  help="Columns to filter on, comma separated")
	opt_parser.add_option("-y", "--filter_cols_vals",
						  dest="filter_cols_vals",
						  help="Values to keep, one per column")
	return opt_parser


def MakeBinary(a):
	return a.lower() == 'true'

def MergeBinary(vals, names):
	if not vals:
		logging.error('Empty vals iterable')
		return None
	
	if len(vals) == 1:
		if MakeBinary(vals[0]):
			return names[0]
		return 'NONE'
	
	bools = map(MakeBinary, vals)
	all_true = reduce(lambda x,y: x and y, bools)
	if all_true:
		if len(vals) == 2:
			return 'BOTH'
		return 'ALL'
	
	true_names = [names[i] for i, test in enumerate(bools)
				  if test == True]
	if true_names:
		return ', '.join(true_names)

	if len(vals) == 2:
		return 'NEITHER'
	return 'NONE'


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

	# Set up filters	
	filter_cols, filter_cols_vals = [], []
	if options.filter_cols:
		assert options.filter_cols_vals
		filter_cols = map(str.strip, options.filter_cols.split(','))
		filter_cols_vals = map(str.strip, options.filter_cols_vals.split(','))

	# Grab relevant data
	pairmap = {}
	all_a = set()
	all_b = set()
	all_vals = set()
	for row in r:
		a, b = row[first_col].strip(), row[second_col].strip()
		if not a or not b:
			continue

		# Apply column filters
		apply_filter = lambda x, y: x in row and row[x] == y
		applied_filters = map(apply_filter, filter_cols, filter_cols_vals)
		or_reduce = lambda x, y: x or y
		passed_filter = reduce(or_reduce, applied_filters, False)
		if filter_cols and not passed_filter:
			continue
		
		all_a.add(a)
		all_b.add(b)

		vals = [b]
		names = [second_col_key]		
		if third_col:
			vals.append(row[third_col].strip())
			names.append(third_col_key)
			
		# Apply value filters
		value = MergeBinary(vals, names)
		if filter_values and value in filter_values:
			continue

		all_vals.add(value)
		pairmap.setdefault(a, []).append(value)
	
	# Compute counts an weights
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

	# Plot bar chart
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

	# Title, ticks, labels
	title = options.title or '%s vs. %s' % (options.first_col, options.second_col)
	pylab.title(title)
	pylab.xlabel(first_col, fontsize='large')

	size_12 = FontProperties(size=12)
	a_labels = [a or "None given" for a in all_a]
	pylab.xticks(ind + 0.25, a_labels, fontproperties=size_12)
	axes.yaxis.set_major_locator(pylab.NullLocator())

	# Show figure.
	pylab.legend()
	pylab.show()
	
	if options.output_filename:
		fieldnames = [''] + a_labels
		f = open(options.output_filename, 'w')
		w = csv.writer(f)
		w.writerow(fieldnames)
		for key in sorted(all_vals):
			row = [key] + [c.get(key, 0) for c in count_array]
			w.writerow(row)
		f.close()


if __name__ == '__main__':
	Main()
