
import json
import sys

from optparse import OptionParser
from pygibbs.kegg import Kegg
from toolbox.database import SqliteDatabase
from pygibbs.thermodynamics import PsuedoisomerTableThermodynamics
from pygibbs.nist_verify import LoadAllEstimators


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


def ExportJSONFiles():
    options, _ = MakeOpts().parse_args(sys.argv)
    print "Using the database file: " + options.public_db
    print "Using the NIST table: " + options.nist_table
    print "Saving the data to the CSV file: " + options.output_csv

    db = SqliteDatabase(options.public_db)
    db.Table2CSV(options.output_csv, options.nist_table, write_titles=True)

if __name__ == '__main__':
    ExportJSONFiles()