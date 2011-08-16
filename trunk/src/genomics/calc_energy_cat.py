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


def IsOrganic(energy_tag):
	return ('hetero' in energy_tag or
			'organo' in energy_tag)


def IsInorganic(energy_tag):
	return ('auto' in energy_tag or
			'litho' in energy_tag or
			'photo' in energy_tag or
			'chemosynth' in energy_tag)

def IsC1(energy_tag):
	return ('methan' in energy_tag or
			'methyl' in energy_tag)

def IsMix(energy_tag):
	return ('mixo' in energy_tag)

def GetEnergyCategory(energy_tags):
	org = False
	inorg = False
	mix = False
	c1 = False
	for tag in energy_tags:
		org |= IsOrganic(tag)
		inorg |= IsInorganic(tag)
		mix |= IsMix(tag)
		c1 |= IsC1(tag)

	if mix:
		return 'Both'
	if org and inorg:
		return 'Both'
	if c1 and inorg:
		return 'Both'
	if c1:
		return 'C1'
	if org:
		return 'Organic'
	if inorg:
		return 'Inorganic'
	return ''
	

def Main():
	options, _ = MakeOpts().parse_args(sys.argv)
	assert options.input_filename
	assert options.col

	l = []
	f = open(options.input_filename)
	for row in csv.DictReader(f):
		if options.col in row:
			l.append(row[options.col].lower().split(", "))

	categories = map(GetEnergyCategory, l)
	for cat in categories:
		print cat	


if __name__ == '__main__':
	Main()
