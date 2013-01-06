from argparse import ArgumentParser
from toolbox.tecan import WriteToDatabase, CollectData, CollectDataFromSingleFile,\
    GetExpDate, GetCurrentExperimentID
import sys, os
from toolbox.database import MySQLDatabase, SqliteDatabase
import time

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

def MakeOpts():
    """Returns an OptionParser object with all the default options."""
    parser = ArgumentParser()
    parser.add_argument("-o", "--host",
                        dest="host",
                        default="hldbv02",
                        help="The hostname for the MySQL database")
    parser.add_argument("-t", "--tar_filename",
                        default=None,
                        help="The filename for the .tar.gz results bundle")
    parser.add_argument("-e", "--exp_id",
                        default=None,
                        help="The ID for the imported experiment data")
    parser.add_argument("-x", "--xml_filename",
                        default=None,
                        help="The filename for a single XML result file")
    parser.add_argument("-a", "--xml_dir",
                        default=None,
                        help="The directory from which to import the latest XML results file")
    parser.add_argument("-p", "--plate_id",
                        default=None, type=int,
                        help="The ID of the read plate (necessary when using -x)")
    parser.add_argument("-g", "--generate_exp_id",
                        default=False, action="store_true", 
                        help="Generate a new Experiment ID and add it to the database")
    parser.add_argument('-d', '--debug', action='store_true',
                        default=False, help='debug mode, store results in dummy DB')
    return parser

def GetTimeString():
    return time.strftime('%Y-%m-%d %H:%M:%S', time.localtime())

def get_latest_file(path):
    filelist = os.listdir(path)
    filelist = filter(lambda x: not os.path.isdir(path + x), filelist)
    return max(filelist, key=lambda x: os.stat(path + x).st_mtime)

def main():
    options = MakeOpts().parse_args()

    if options.debug:
        db = CreateDummyDB()
    else:
        db = MySQLDatabase(host=options.host, user='ronm', port=3306,
                           passwd='a1a1a1', db='tecan')

    if options.tar_filename:
        """
            Imports a list of XML files bundled in a TAR file into the database.
            The Experiment ID is determined by the time-stamp of the first measurement.
            If this Experiment ID 
        """
        
        if not os.path.exists(options.tar_filename):
            print "File not found: " + options.tar_filename
            sys.exit(-1)
        
        print "Importing from file: " + options.tar_filename
        MES = CollectData(options.tar_filename)

        exp_id = options.exp_id or GetExpDate(MES)        
        print "Experiment ID: " + exp_id
        
        # delete any previous data regarding this exp_id
        db.Execute("DELETE FROM tecan_readings WHERE exp_id='%s'" % exp_id)
        db.Execute("DELETE FROM tecan_experiments WHERE exp_id='%s'" % exp_id)
        db.Insert('tecan_experiments', [exp_id, options.plate_id, 
                                        "Imported from TAR file on " + 
                                        GetTimeString()])
        WriteToDatabase(MES, db, exp_id)
    elif options.xml_dir or options.xml_filename:
        if options.plate_id is None:
            print "Plate ID not supplied, but is mandatory when importing an XML"
            sys.exit(-1)
    if options.xml_dir:
        if not os.path.exists(options.xml_dir):
            print "Directory not found: " + options.xml_dir
            sys.exit(-1)
        options.xml_filename = options.xml_dir + get_latest_file(options.xml_dir)
        if not os.path.exists(options.xml_filename):
            print "File not found: " + options.xml_filename
            sys.exit(-1)
        
        print "Importing from file: " + options.xml_filename
        MES = CollectDataFromSingleFile(options.xml_filename, options.plate_id)
        print MES
        exp_id = GetCurrentExperimentID(db)
        print "Experiment ID: " + exp_id
        
        WriteToDatabase(MES, db, exp_id)
   
    db.Commit()
    print "Done!"
    
if __name__ == '__main__':
    main()
