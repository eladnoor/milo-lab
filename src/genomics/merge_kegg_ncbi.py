#!/usr/bin/python

import csv
import logging
import sys

from optparse import OptionParser

def MakeOpts():
    """Returns an OptionParser object with all the default options."""
    opt_parser = OptionParser()
    opt_parser.add_option("-k", "--kegg_filename",
                          dest="kegg_filename",
                          help="input CSV with mapping of KEGG species ids to NCBI taxa")
    opt_parser.add_option("-n", "--ncbi_filename",
                          dest="ncbi_filename",
                          help="input CSV with NCBI data")
    opt_parser.add_option("-o", "--output_filename",
                          dest="output_filename",
                          help="filename to write output csv to")
    return opt_parser


def Main():
    options, _ = MakeOpts().parse_args(sys.argv)
    assert options.kegg_filename
    assert options.ncbi_filename
    assert options.output_filename
    print 'Reading KEGG names from', options.kegg_filename
    print 'Reading NCBI data from', options.ncbi_filename

    rows_by_taxa = {}
    r = csv.DictReader(open(options.ncbi_filename), dialect=csv.excel_tab)
    ncbi_fieldnames = list(r.fieldnames)
    for row in r:
        taxa = row.get('NCBI Taxon ID', None)
        if not taxa:
            logging.warning('Undefined taxa for row %s', row)
            continue
        if taxa in rows_by_taxa:
            logging.warning('Duplicate taxa %s', taxa)
            continue

        try:
            stripped = taxa.strip()
            num_taxa = int(taxa.strip())
            rows_by_taxa[stripped] = row
        except Exception, e:
            print e
            print row
            continue
        
    kegg_by_taxa = {}
    for row in csv.DictReader(open(options.kegg_filename)):
        taxa = row.get('NCBI Taxon ID', None)
        kegg_id = row.get('KEGG Organism ID')
        if taxa is not None:
            if taxa not in kegg_by_taxa:
                stripped = taxa.strip()
                kegg_by_taxa[stripped] = kegg_id
            else:
                logging.warning('Duplicate taxa %s', taxa)
                

    print 'Writing output to', options.output_filename
    ncbi_fieldnames.append('KEGG ID')
    w = csv.DictWriter(open(options.output_filename, 'w'), ncbi_fieldnames)
    w.writeheader()

    for taxa, row in rows_by_taxa.iteritems():
        kegg_id = kegg_by_taxa.get(taxa, None)
        
        row['KEGG ID'] = kegg_id
        w.writerow(row)


if __name__ == '__main__':
	Main()
