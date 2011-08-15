#!/usr/bin/python

import csv
import itertools
import pylab
import sys

from matplotlib.font_manager import FontProperties
from optparse import OptionParser


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
	opt_parser.add_option("-f", "--filter_col",
						  dest="filter_col",
						  help="A column to filter based on.")
	opt_parser.add_option("-v", "--filter_val",
						  dest="filter_val",
						  help="The value to filter out.")
	return opt_parser


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

	# Read and filter species data
	r = csv.DictReader(open(options.species_filename))
	filter_col, filter_val = options.filter_col, options.filter_val
	first_col, second_col = options.first_col, options.second_col
	pairmap = {}
	for row in r:
		if (filter_col and filter_col in row and
			row[filter_col] == filter_val):
			continue

		key = (row[first_col], row[second_col])
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

	# Count instances of pairs.
	counts = dict((k, len(v)) for k,v in pairmap.iteritems())

	# Format for pylab plotting.
	x_vals = []
	y_vals = []
	z_vals = []
	for pair in all_possible_pairs:
		a,b = pair
		x_vals.append(a_to_num[a])
		y_vals.append(b_to_num[b])
		z_vals.append(counts.get(pair, 0))
	
	# Scale counts to radii, max 0.5 so they don't overlap with adjacent blobs.
	scaled_z_vals = pylab.array(map(float, z_vals))
	max_z = max(scaled_z_vals)
	scaled_z_vals /= (2.0*max_z)
	
	# Plot circle scatter.
	axes = pylab.axes()
	axes.grid(color='g', linestyle='--', linewidth=1)
	CircleScatter(axes, x_vals, y_vals, scaled_z_vals)

	# Labels, titles and ticks.
	pylab.title('%s vs. %s' % (options.first_col, options.second_col))
	pylab.xlabel(options.first_col)
	pylab.ylabel(options.second_col)
	pylab.figtext(0.70, 0.02, '%d examples total.' % sum(counts.values()))

	size_8 = FontProperties(size=8)
	a_labels = [a or "None given" for a in all_a]
	b_labels = [b or "None given" for b in all_b]
	pylab.xticks(range(0, len(a_labels)), a_labels,
				 fontproperties=size_8)
	pylab.yticks(range(0, len(b_labels)), b_labels,
				 fontproperties=size_8)
	for x, y, z in zip(x_vals, y_vals, z_vals):
		pylab.text(x+0.1, y+0.1, str(z))

	# Scale and show.
	pylab.axis('scaled')
	pylab.axis([-1,len(all_a),-1,len(all_b)])
	pylab.show()
	

if __name__ == '__main__':
	Main()

