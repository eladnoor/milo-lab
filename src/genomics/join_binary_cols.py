#!/usr/bin/python

import csv
import logging
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
	opt_parser.add_option("-a", "--col_a",
						  dest="col_a",
						  help="The first column.")
	opt_parser.add_option("-b", "--col_b",
						  dest="col_b",
						  help="The second column.")
	opt_parser.add_option("-c", "--col_c",
						  dest="col_c",
						  help="The third, output column.")
	opt_parser.add_option("-v", "--val_a",
						  dest="val_a",
						  help="The value to use when col a is true only.")
	opt_parser.add_option("-w", "--val_b",
						  dest="val_b",
						  help="The value to use when col b is true only.")
	return opt_parser


def ParseBool(v):
	l = v.lower()
	return l == 'true'


def Main():
	options, _ = MakeOpts().parse_args(sys.argv)
	assert options.input_filename and options.output_filename
	assert options.col_a and options.col_b and options.col_c

	val_a, val_b = options.val_a, options.val_b
	col_a, col_b, col_c = options.col_a, options.col_b, options.col_c
	r = csv.DictReader(open(options.input_filename))
	fieldnames = r.fieldnames
	fieldnames.append(col_c)
	w = csv.DictWriter(open(options.output_filename, 'w'), fieldnames)
	w.writeheader()
	for row in r:
		
		if (col_a not in row or
			col_b not in row or
			col_c not in row):
			logging.warning('Missing data for col.')
			continue

		a = ParseBool(row[col_a])
		b = ParseBool(row[col_b])
		if a and b:
			row[col_c] = 'BOTH'
		elif a:
			row[col_c] = val_a
		elif b:
			row[col_c] = val_b
			
		w.writerow(row)


if __name__ == '__main__':
	Main()
