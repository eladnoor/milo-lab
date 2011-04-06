#!/usr/bin/python

import sys

from optparse import OptionParser
from pygibbs.thermodynamics import PsuedoisomerTableThermodynamics
from pygibbs.kegg import Kegg
from pygibbs.groups import GroupContribution
from toolbox.database import SqliteDatabase


def MakeOpts():
    """Returns an OptionParser object with all the default options."""
    opt_parser = OptionParser()
    opt_parser.add_option("-p", "--public_database_location", dest="public_db_filename",
                          default="../data/public_data.sqlite",
                          help="The public database location")
    opt_parser.add_option("-d", "--database_location", dest="db_filename",
                          default="../res/gibbs.sqlite",
                          help="The database location")
    opt_parser.add_option("-o", "--out_filename",
                          dest="out_filename",
                          default='../res/pseudoisomers.json',
                          help="Output filename.")
    opt_parser.add_option("-t", "--thermodynamics_filename",
                          dest="thermo_filename",
                          default='../data/thermodynamics/dG0.csv',
                          help="Thermodynamics filename.")
    
    return opt_parser


def main():
    # References:
    # Alberty 2006: Alberty R. A. - Biochemical Thermodynamics: Applications of Mathematica (Methods of Biochemical Analysis), Wiley 2006
    # 
    options, _ = MakeOpts().parse_args(sys.argv)
    print 'Database filename:', options.db_filename
    print 'Thermodynamics filename:', options.thermo_filename
    print 'Output JSON filename:', options.out_filename

    db = SqliteDatabase(options.db_filename)
    db_public = SqliteDatabase(options.public_db_filename)
    
    alberty = PsuedoisomerTableThermodynamics.FromDatabase(db_public, 'alberty_pseudoisomers')
    group_contrib = PsuedoisomerTableThermodynamics.FromDatabase(db, 'gc_pseudoisomers')
    group_contrib.override_data(alberty)
    group_contrib.write_data_to_json(options.out_filename)
    print 'Done.'

    
if __name__ == "__main__":
    main()
