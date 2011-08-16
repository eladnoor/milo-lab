#!/usr/bin/python

import csv
import sys

from optparse import OptionParser


def MakeOpts():
	"""Returns an OptionParser object with all the default options."""
	opt_parser = OptionParser()
	opt_parser.add_option("-i", "--input_filename",
						  dest="input_filename",
						  help="input CSV")
	opt_parser.add_option("-o", "--output_filename",
						  dest="output_filename",
						  help="output CSV")
	opt_parser.add_option("-f", "--filter_cols",
						  dest="filter_cols",
						  help="The columns to filter based on (comma separated)")
	opt_parser.add_option("-v", "--filter_vals",
						  dest="filter_vals",
						  help="The values to of examples to keep (comma separated).")
	return opt_parser


def Main():
	options, _ = MakeOpts().parse_args(sys.argv)
	assert options.input_filename and options.output_filename
	assert options.filter_cols and options.filter_vals

	filter_cols = map(str.strip, options.filter_cols.split(','))
	filter_vals = map(str.strip, options.filter_vals.split(','))
	assert len(filter_cols) == len(filter_vals)

	r = csv.DictReader(open(options.input_filename))
	rows = []
	
	for row in r:
		apply_filter = lambda x, y: x in row and row[x] == y
		or_reduce = lambda x, y: x or y
		passed_filter = reduce(or_reduce, map(apply_filter,
											  filter_cols, filter_vals))
		if passed_filter:
			rows.append(row)


	w = csv.DictWriter(open(options.output_filename, 'w'), r.fieldnames)
	w.writeheader()
	map(w.writerow, rows)

	
if __name__ == '__main__':
	Main()

