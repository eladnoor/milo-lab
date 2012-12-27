from argparse import ArgumentParser
import sys, os
from toolbox.database import MySQLDatabase, SqliteDatabase
import time
from pyvictor.victor_parser import VictorParser

# all Victor experiments have only one plate, so their ID is 0
PLATE_ID = 0

def MakeOpts():
    """Returns an OptionParser object with all the default options."""
    parser = ArgumentParser()
    parser.add_argument('-d', '--debug', action='store_true',
                        help='debug mode, store results in dummy DB')
    parser.add_argument('-e', '--exp_id', default=None,
                        help='Override the experiment ID '
                             '(default is to use the measurement date-time)')
    parser.add_argument('xls_file',
                        help='The path to the XLS file containing Victor results')
    return parser

def GetTimeString():
    return time.strftime('%Y-%m-%d %H:%M:%S', time.localtime())

def GetLatestFile(path):
    filelist = [f for f in os.listdir(path) if f[-4:] == '.xls']
    filelist = filter(lambda x: not os.path.isdir(path + x), filelist)
    return max(filelist, key=lambda x: os.stat(path + x).st_mtime)

def CreateDummyDB():
    db = SqliteDatabase('/tmp/victor_dummy.sqlite', 'w')
    db.CreateTable('tecan_readings', 'exp_id TEXT, plate_id TEXT, '
                   'measurement_name TEXT, row INT, col INT, time_in_sec INT, '
                   'value REAL', drop_if_exists=False)
    db.CreateTable('tecan_experiments', 'exp_id TEXT, desc TEXT',
                   drop_if_exists=False)
    db.CreateTable('tecan_plates', 'exp_id TEXT, plate_id TEXT, '
                   'field1 TEXT, field2 INT, field3 INT', drop_if_exists=False)
    return db

def main():
    """
        Imports an XLS file of Victor exported results.
        The Experiment ID is determined by the time-stamp of the first measurement.
    """
    
    options = MakeOpts().parse_args()

    if options.debug:
        db = CreateDummyDB()
    else:
        db = MySQLDatabase(host='hldbv02', user='ronm', port=3306,
                           passwd='a1a1a1', db='tecan')

    if options.xls_file is not None:
        xls_filename = options.xls_file
    else:
        if not os.path.exists(options.xls_dir):
            print "Directory not found: " + options.xls_dir
            sys.exit(-1)
        xls_filename = options.xls_dir + GetLatestFile(options.xls_dir)
        
    if not os.path.exists(xls_filename):
        print "File not found: " + xls_filename
        sys.exit(-1)
    
    print "Importing from file: " + xls_filename
    vp = VictorParser()
    vp.parse_excel(open(xls_filename, 'r'))
    
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
    db.Insert('tecan_plates', [exp_id, PLATE_ID, "", None, None])
    vp.write_to_database(db, exp_id)
    db.Commit()
    
    if options.debug:
        print "Done, go check out the results at %s" % db.filename
    else:
        print "Done, go check out the results at http://eladpc1/RoboSite/Exp/%s/0" % exp_id
    
if __name__ == '__main__':
    main()
