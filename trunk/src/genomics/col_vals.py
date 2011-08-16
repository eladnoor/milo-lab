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
	opt_parser.add_option("-c", "--col",
						  dest="col",
						  help="column to grab")
	return opt_parser

def Main():
	options, _ = MakeOpts().parse_args(sys.argv)
	assert options.input_filename
	assert options.col


	s = set()
	f = open(options.input_filename)
	for row in csv.DictReader(f):
		if options.col in row:
			s.update(row[options.col].split(", "))

	print '%d values' % len(s)
	for val in sorted(s):
		print val

if __name__ == '__main__':
	Main()

	
