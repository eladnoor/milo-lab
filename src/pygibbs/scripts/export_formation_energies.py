#!/usr/bin/python

import sys

from optparse import OptionParser
from pygibbs.thermodynamics import PsuedoisomerTableThermodynamics
from toolbox.database import SqliteDatabase

def MakeOpts():
    """Returns an OptionParser object with all the default options."""
    opt_parser = OptionParser()
    opt_parser.add_option("-d", "--database_location", dest="db_filename",
                          default="../res/gibbs.sqlite",
                          help="The database location")
    opt_parser.add_option("-o", "--out_filename",
                          dest="out_filename",
                          default='../res/pseudoisomers.json',
                          help="Output filename")
    opt_parser.add_option("-t", "--thermodynamics_filename",
                          dest="thermo_filename",
                          default='../data/thermodynamics/dG0.csv',
                          help="Thermodynamics filename")
    opt_parser.add_option("-g", "--gc_table_name",
                          dest="gc_table_name",
                          default='gc_pseudoisomers',
                          help="Group Contribution Table Name")
    
    return opt_parser


def main():
    # References:
    # Alberty 2006: Alberty R. A. - Biochemical Thermodynamics: Applications of Mathematica (Methods of Biochemical Analysis), Wiley 2006
    # 
    options, _ = MakeOpts().parse_args(sys.argv)
    print 'Database filename:', options.db_filename
    print 'Thermodynamics filename:', options.thermo_filename
    print 'Output JSON filename:', options.out_filename
    print 'Group Contribution Table Name:', options.gc_table_name

    db = SqliteDatabase(options.db_filename)
    
    observed_thermo = PsuedoisomerTableThermodynamics.FromCsvFile(
                                                options.thermo_filename)
    if not db.DoesTableExist(options.gc_table_name):
        raise ValueError('The table %s does not exist in the database. '
                         'Please run the groups.py script and try again.'
                         % options.gc_table_name)
        
    group_contrib = PsuedoisomerTableThermodynamics.FromDatabase(
                                                db, options.gc_table_name)
    group_contrib.override_data(observed_thermo)
    group_contrib.write_data_to_json(options.out_filename)
    print 'Done.'

    
if __name__ == "__main__":
    main()
