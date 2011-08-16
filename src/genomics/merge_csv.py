#!/usr/bin/python

import csv
import logging
import sys

from optparse import OptionParser


def MakeOpts():
	"""Returns an OptionParser object with all the default options."""
	opt_parser = OptionParser()
	opt_parser.add_option("-a", "--file_a",
						  dest="file_a",
						  help="First CSV file.")
	opt_parser.add_option("-b", "--file_b",
						  dest="file_b",
						  help="Second CSV file.")
	opt_parser.add_option("-k", "--merge_key",
						  dest="merge_key",
						  help="The column to merge on.")	
	opt_parser.add_option("-o", "--output_file",
						  dest="output_file",
						  help="Output CSV filename.")
	return opt_parser


def CsvToDict(fname, key_col):
	d = {}
	f = open(fname)
	r = csv.DictReader(f)
	for row in r:
		if key_col not in row:
			logging.warning('Row missing value for "%s"', key_col)
			continue

		key = row.get(key_col, None)
		if not key:
			logging.debug('Skipping row without value for "%s"', key_col)
			continue

		if key in d:
			logging.warning('Duplicate %s: %s', key_col, key)
			continue

		d[key] = row

	f.close()
	return d, r.fieldnames


def MergeCsvDicts(da, db):
	"""Merge two csv files, defaulting to the first."""
	overlap = set(da.keys()).intersection(db.keys())
	for key in overlap:
		d = dict(da[key])
		d.update(db[key])
		yield d


def Main():
	options, _ = MakeOpts().parse_args(sys.argv)
	assert options.file_a and options.file_b
	assert options.merge_key
	assert options.output_file
	print 'Reading files', options.file_a, 'and', options.file_b
	print 'Merging on column', options.merge_key

	da, fieldsa = CsvToDict(options.file_a, options.merge_key)
	db, fieldsb = CsvToDict(options.file_b, options.merge_key)

	fieldnames = sorted(set(fieldsa + fieldsb))
	outfile = open(options.output_file, 'w')
	w = csv.DictWriter(outfile, fieldnames)
	w.writeheader()
	map(w.writerow, MergeCsvDicts(da, db))


if __name__ == '__main__':
	Main()


