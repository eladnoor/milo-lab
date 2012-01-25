from optparse import OptionParser
import sys, os
from toolbox.database import MySQLDatabase
import time
from pyvictor.victor_parser import VictorParser

def MakeOpts():
    """Returns an OptionParser object with all the default options."""
    opt_parser = OptionParser(usage="usage: %prog [options] xls_fname")
    opt_parser.add_option("-o", "--host",
                          dest="host",
                          default="hldbv02",
                          help="The hostname for the MySQL database (default=%default)")
    opt_parser.add_option("-e", "--exp_id",
                          dest="exp_id",
                          default=None,
                          help="The ID for the imported experiment data, "
                          "default is to use the measurement date.")
    opt_parser.add_option("-p", "--plate_id",
                          dest="plate_id", default=0, type='int',
                          help="The ID of the read plate (should always be 0)")
    return opt_parser

def GetTimeString():
    return time.strftime('%Y-%m-%d %H:%M:%S', time.localtime())

def main():
    """
        Imports an XLS file of Victor exported results.
        The Experiment ID is determined by the time-stamp of the first measurement.
    """
    
    opt_parser = MakeOpts()
    options, args = opt_parser.parse_args()

    db = MySQLDatabase(host=options.host, user='ronm', port=3306,
                       passwd='a1a1a1', db='tecan')

    if not args:
        print opt_parser.get_usage()
        sys.exit(-1)
    xls_filename = args[0]
    if not os.path.exists(xls_filename):
        print "File not found: " + xls_filename
        sys.exit(-1)
    
    print "Importing from file: " + xls_filename
    vp = VictorParser()
    vp.parse_excel(xls_filename)
    
    exp_id = options.exp_id or vp.get_time_string()   
    print "Experiment ID: " + exp_id
    
    # delete any previous data regarding this exp_id
    db.Execute("DELETE FROM tecan_readings WHERE exp_id='%s'" % exp_id)
    db.Execute("DELETE FROM tecan_experiments WHERE exp_id='%s'" % exp_id)
    db.Insert('tecan_experiments', [exp_id, "Imported from XLS file on " + 
                                    GetTimeString()])
    db.Insert('tecan_plates', [exp_id, options.plate_id, ""])
    vp.write_to_database(db, exp_id)
    db.Commit()
    print "Done!"
    
if __name__ == '__main__':
    main()
