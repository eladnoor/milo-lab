from optparse import OptionParser
from pytecan.util import WriteToSqlite, CollectData
import sys
#import sqlite3
import MySQLdb as mdb
from toolbox.database import MySQLDatabase


def MakeOpts():
    """Returns an OptionParser object with all the default options."""
    opt_parser = OptionParser()
    #opt_parser.add_option("-d", "--sqlite_db_filename",
    #                      dest="sqlite_db_filename",
    #                      default="../res/tecan.sqlite",
    #                      help="The filename of the Sqlite database")
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
    #print "Importing into database: " + options.sqlite_db_filename
    if options.exp_id:
        print "Experiment ID: " + options.exp_id
    
    MES = CollectData(options.tar_filename)
    db = MySQLDatabase(host='eladpc1', user='eladn', passwd='a1a1a1', db='tecan')
    exp_id = WriteToSqlite(MES, db, options.exp_id)
    db.Commit()
    print "Done importing experiment: %s" % exp_id
    
    for row in db.Execute('SELECT exp_id, plate, count(*) from tecan_readings '
                          'GROUP BY exp_id, plate'):
        print row
    
if __name__ == '__main__':
    main()