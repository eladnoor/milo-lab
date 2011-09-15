from optparse import OptionParser
import sys
from toolbox.database import MySQLDatabase
import time

def MakeOpts():
    """Returns an OptionParser object with all the default options."""
    opt_parser = OptionParser()
    #opt_parser.add_option("-d", "--sqlite_db_filename",
    #                      dest="sqlite_db_filename",
    #                      default="../res/tecan.sqlite",
    #                      help="The filename of the Sqlite database")
    opt_parser.add_option("-p", "--plate_num",
                          type='int',
                          dest="plate_num",
                          default=1,
                          help="The number for the plate that is to be exported")
    opt_parser.add_option("-e", "--exp_id",
                          dest="exp_id",
                          default=None,
                          help="The expID for the data")
    return opt_parser

def main():
    opt_parser = MakeOpts()
    options, _ = opt_parser.parse_args(sys.argv)
    if not options.exp_id:
        opt_parser.print_help(sys.stderr)
        sys.exit(-1)

    #print "Importing into database: " + options.sqlite_db_filename
    print "Experiment ID: " + options.exp_id
    print "Plate num: %d" % options.plate_num
    
    db = MySQLDatabase(host='132.77.80.238', user='ronm', 
                       passwd='a1a1a1', db='tecan')
    
    for row in db.Execute("SELECT * from tecan_readings "
                          "WHERE exp_id='%s' AND plate='%d' "
                          "ORDER BY exp_id, plate, reading_label, row, col;" %
                          (options.exp_id, options.plate_num)):
        (_exp_id, _plate, reading_label, row, col, time_in_sec, measurement) = row
        t = time.strftime('%Y-%m-%d %H:%M:%S', time.gmtime(time_in_sec))
        print "(%d,%d) %s, %s : %.4f" % (row, col, reading_label, t, measurement)
    
    del db
    
if __name__ == '__main__':
    main()