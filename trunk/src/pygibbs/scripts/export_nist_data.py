
import json
import sys

from optparse import OptionParser
from pygibbs.kegg import Kegg
from toolbox.database import SqliteDatabase
from pygibbs.thermodynamics import PsuedoisomerTableThermodynamics
from pygibbs.nist_verify import LoadAllEstimators
import csv


def MakeOpts():
    """Returns an OptionParser object with all the default options."""
    opt_parser = OptionParser()
    opt_parser.add_option("-n", "--nist_table",
                          dest="nist_table",
                          default="nist_equilibrium",
                          help="the name of the NIST table")
    opt_parser.add_option("-d", "--public_db",
                          dest="public_db",
                          default="../data/public_data.sqlite",
                          help="the name of the DB file that contains the NIST data")
    opt_parser.add_option("-o", "--output_csv",
                          dest="output_csv",
                          default='../res/nist_equilibrium.csv',
                          help="the name of the output CSV file")
    return opt_parser


def reformat_number_string(s, new_format):
    if s:
        return new_format % float(s)
    else:
        return None

def ExportJSONFiles():
    options, _ = MakeOpts().parse_args(sys.argv)
    print "Using the database file: " + options.public_db
    print "Using the NIST table: " + options.nist_table
    print "Saving the data to the CSV file: " + options.output_csv

    db = SqliteDatabase(options.public_db)
    csv_writer = csv.writer(open(options.output_csv, 'w'))
    csv_writer.writerow(['url','reference_id','method','evaluation','ec',
                         'enzyme','kegg_reaction','reaction','K','K_tag',
                         'T (K)','I (M)','pH','pMg'])
    for row in db.DictReader(options.nist_table):
        csvrow = [row[t] for t in ['url','reference_id','method','evaluation','ec','enzyme','kegg_reaction','reaction']]
        csvrow += [reformat_number_string(row['K'], '%.3e'), 
                   reformat_number_string(row['K_tag'], '%.3e'),
                   reformat_number_string(row['T'], '%.2f'),
                   reformat_number_string(row['I'], '%.2f'),
                   reformat_number_string(row['pH'], '%.2f'),
                   reformat_number_string(row['pMg'], '%.2f')]
        csv_writer.writerow(csvrow)

if __name__ == '__main__':
    ExportJSONFiles()