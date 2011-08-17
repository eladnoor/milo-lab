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
	opt_parser.add_option("-s", "--species_filename",
						  dest="species_filename",
						  help="input CSV with species names and IDs")
	opt_parser.add_option("-a", "--first_col",
						  dest="first_col",
						  help="First column in 2-way histogram")
	opt_parser.add_option("-b", "--second_col",
						  dest="second_col",
						  help="First column in 2-way histogram")
	opt_parser.add_option("-c", "--third_col",
						  dest="third_col",
						  help="THIRD DIMENSION!")
	opt_parser.add_option("-f", "--filter_cols",
						  dest="filter_cols",
						  help="The columns to filter based on (comma separated)")
	opt_parser.add_option("-v", "--filter_vals",
						  dest="filter_vals",
						  help="The values to of examples to keep (comma separated).")
	opt_parser.add_option("-l", "--limit_val",
						  dest="limit_val",
						  help="limit analysis to examples where at least one column has this value.")
	return opt_parser


COLORS = 'grbcmykw'


def PieScatter(axes, x_array, y_array, counts_array, max_count, colormap, **kwargs):
    for x, y, counts in zip(x_array, y_array, counts_array):
		total_count = sum(counts.itervalues())
		radius = pylab.sqrt(float(total_count) / (5.0 * float(max_count)))
		axes.text(x + 0.05, y - radius - 0.1, str(total_count))

		theta_1 = 0
		for i, (key, value) in enumerate(counts.iteritems()):
			color = colormap[key]
			sweep = (float(value) / float(total_count)) * 360.0
			w = Wedge((x,y), radius, theta_1, theta_1 + sweep, color=color, label=key)
			theta_1 += sweep
			axes.add_patch(w)

    return True


def CircleScatter(axes, x_array, y_array, z_array, **kwargs):
    for x, y, radius in zip(x_array, y_array, z_array):
        circle = pylab.Circle((x,y), radius=radius, **kwargs)
        axes.add_patch(circle)
    return True


def Main():
	options, _ = MakeOpts().parse_args(sys.argv)
	assert options.species_filename
	assert options.first_col and options.second_col
	print 'Reading species list from', options.species_filename

	filter_cols, filter_vals = [], []
	if options.filter_cols:
		assert options.filter_vals
		filter_cols = map(str.strip, options.filter_cols.split(','))
		filter_vals = map(str.strip, options.filter_vals.split(','))
	

	# Read and filter species data
	r = csv.DictReader(open(options.species_filename))
	first_col, second_col = options.first_col, options.second_col
	pairmap = {}
	for row in r:
		apply_filter = lambda x, y: x in row and row[x] == y
		or_reduce = lambda x, y: x or y
		passed_filter = reduce(or_reduce, map(apply_filter,
											  filter_cols, filter_vals), True)
		if not passed_filter:
			continue
		
		a, b = row[first_col].strip(), row[second_col].strip()
		if not a or not b:
			continue
		key = (a, b)
		pairmap.setdefault(key, []).append(row)

	#for key, row_list in pairmap.iteritems():
	#	if len(row_list) < 5:
	#		print key, row_list

	# Find cross-product of column values (all pairs).
	all_a = list(set([x[0] for x in pairmap.keys()]))
	all_b = list(set([x[1] for x in pairmap.keys()]))
	a_to_num = dict((v,i) for i,v in enumerate(all_a))
	b_to_num = dict((v,i) for i,v in enumerate(all_b))
	all_possible_pairs = list(itertools.product(all_a, all_b))

	third_col = options.third_col or 'fake key'
	get_col = lambda x: x.get(third_col, None)
	col_vals = dict((k, map(get_col, v)) for k, v in pairmap.iteritems())
	counts = {}
	totals = []
	all_vals = set()
	for k, v in col_vals.iteritems():
		counter = Counter(v)
		all_vals.update(counter.keys())
		counts[k] = counter
		totals.append(sum(counter.values()))

	x_vals = []
	y_vals = []
	count_array = []
	z_vals = []
	max_val = max(totals)
	for pair in all_possible_pairs:
		a,b = pair
		x_vals.append(a_to_num[a])
		y_vals.append(b_to_num[b])
		z_vals.append(sum(counts.get(pair, {}).values()))
		count_array.append(counts.get(pair, {}))
		
	# Plot circle scatter.
	axes = pylab.axes()
	axes.grid(color='g', linestyle='--', linewidth=1)
	
	if options.third_col:
		colormap = ColorMap(all_vals)
		PieScatter(axes, x_vals, y_vals, count_array, max_val, colormap)

		handles, labels = axes.get_legend_handles_labels()
		mapped_labels = dict(zip(labels, handles))
		labels = sorted(mapped_labels.keys())
		handles = [mapped_labels[k] for k in labels]
		pylab.legend(handles, labels)
	else:
		scaled_z_vals = pylab.array(map(float, z_vals))
		max_z = max(scaled_z_vals)
		scaled_z_vals /= (4.0*max_z)
		scaled_z_vals = pylab.sqrt(scaled_z_vals)

		CircleScatter(axes, x_vals, y_vals, scaled_z_vals)
		for x, y, z in zip(x_vals, y_vals, z_vals):
			pylab.text(x+0.1, y+0.1, str(z))
		

	# Labels, titles and ticks.
	pylab.title('%s vs. %s' % (options.first_col, options.second_col))
	pylab.xlabel(options.first_col)
	pylab.ylabel(options.second_col)
	pylab.figtext(0.70, 0.02, '%d examples total.' % sum(totals))

	size_8 = FontProperties(size=8)
	a_labels = [a or "None given" for a in all_a]
	b_labels = [b or "None given" for b in all_b]
	pylab.xticks(range(0, len(a_labels)), a_labels,
				 fontproperties=size_8)
	pylab.yticks(range(0, len(b_labels)), b_labels,
				 fontproperties=size_8)


	# Scale and show.
	pylab.axis('scaled')
	pylab.axis([-1,len(all_a),-1,len(all_b)])
	pylab.show()
	

if __name__ == '__main__':
	Main()

