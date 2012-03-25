from optparse import OptionParser
import sys, os
from toolbox.database import MySQLDatabase
import time
from pyvictor.victor_parser import VictorParser

def MakeOpts():
    """Returns an OptionParser object with all the default options."""
    opt_parser = OptionParser(usage="usage: %prog [options]")
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
    opt_parser.add_option("-d", "--xls_dir",
                          dest="xls_dir", default='/media/vicky_elad/',
                          help="The path to the directory containing Victor XLS results")
    opt_parser.add_option("-i", "--xls_file",
                          dest="xls_file", default=None,
                          help="The path to the XLS file containing Victor results")
    return opt_parser

def GetTimeString():
    return time.strftime('%Y-%m-%d %H:%M:%S', time.localtime())

def get_latest_file(path):
    filelist = [f for f in os.listdir(path) if f[-4:] == '.xls']
    filelist = filter(lambda x: not os.path.isdir(path + x), filelist)
    return max(filelist, key=lambda x: os.stat(path + x).st_mtime)

def main():
    """
        Imports an XLS file of Victor exported results.
        The Experiment ID is determined by the time-stamp of the first measurement.
    """
    
    opt_parser = MakeOpts()
    options, _ = opt_parser.parse_args()

    db = MySQLDatabase(host=options.host, user='ronm', port=3306,
                       passwd='a1a1a1', db='tecan')

    if options.xls_file is not None:
        xls_filename = options.xls_file
    else:
        if not os.path.exists(options.xls_dir):
            print "Directory not found: " + options.xls_dir
            sys.exit(-1)
        xls_filename = options.xls_dir + get_latest_file(options.xls_dir)
        
    if not os.path.exists(xls_filename):
        print "File not found: " + xls_filename
        sys.exit(-1)
    
    print "Importing from file: " + xls_filename
    vp = VictorParser()
    vp.parse_excel(xls_filename)
    
    exp_id = options.exp_id or vp.get_time_string()   
    print "Experiment ID: " + exp_id

    if raw_input('Ready to import Victor results? [y/n] ') != 'y':
        sys.exit(0)
    
    # delete any previous data regarding this exp_id
    db.Execute("DELETE FROM tecan_readings WHERE exp_id='%s'" % exp_id)
    db.Execute("DELETE FROM tecan_experiments WHERE exp_id='%s'" % exp_id)
    db.Execute("DELETE FROM tecan_plates WHERE exp_id='%s'" % exp_id)
    
    desc = "Imported from Victor on " + GetTimeString()
    db.Insert('tecan_experiments', [exp_id, desc])
    db.Insert('tecan_plates', [exp_id, options.plate_id, "", None, None])
    vp.write_to_database(db, exp_id)
    db.Commit()
    print "Done, go check out the results at http://eladpc1/RoboSite/Exp/%s/0" % exp_id
    
if __name__ == '__main__':
    main()
