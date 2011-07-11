from optparse import OptionParser
from pytecan.util import WriteToSqlite, CollectData
import sys
import sqlite3


def MakeOpts():
    """Returns an OptionParser object with all the default options."""
    opt_parser = OptionParser()
    opt_parser.add_option("-d", "--sqlite_db_filename",
                          dest="sqlite_db_filename",
                          default="../res/tecan.sqlite",
                          help="The filename of the Sqlite database")
    opt_parser.add_option("-t", "--tar_filename",
                          dest="tar_filename",
                          default=None,
                          help="The filename for the .tar.gz results bundle")
    opt_parser.add_option("-e", "--exp_id",
                          dest="exp_id",
                          default=None,
                          help="The ID for the imported experiment data")
    return opt_parser

def main():
    opt_parser = MakeOpts()
    options, _ = opt_parser.parse_args(sys.argv)
    if not options.tar_filename:
        opt_parser.print_help(sys.stderr)
        sys.exit(-1)
    print "Importing from file: " + options.tar_filename
    print "Importing into database: " + options.sqlite_db_filename
    if options.exp_id:
        print "Experiment ID: " + options.exp_id
    
    MES = CollectData(options.tar_filename)  
    comm = sqlite3.connect(options.sqlite_db_filename)
    WriteToSqlite(MES, comm, options.exp_id)
    
if __name__ == '__main__':
    main()