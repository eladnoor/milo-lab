from optparse import OptionParser
from pytecan.util import WriteToDatabase, CollectData, CollectDataFromSingleFile,\
    GetExpDate, GetCurrentExperimentID
import sys, os
from toolbox.database import MySQLDatabase
import time

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
    opt_parser.add_option("-x", "--xml_filename",
                          dest="xml_filename",
                          default=None,
                          help="The filename for a single XML result file")
    opt_parser.add_option("-p", "--plate_id",
                          dest="plate_id", default=None, type='int',
                          help="The ID of the read plate (necessary when using -x)")
    opt_parser.add_option("-g", "--generate_exp_id",
                          dest="generate_exp_id", default=False, action="store_true", 
                          help="Generate a new Experiment ID and add it to the database")
    return opt_parser

def GetTimeString():
    return time.strftime('%Y-%m-%d %H:%M:%S', time.localtime())

def main():
    opt_parser = MakeOpts()
    options, _ = opt_parser.parse_args(sys.argv)

    db = MySQLDatabase(host='132.77.80.238', user='ronm', 
                       passwd='a1a1a1', db='tecan')

    if options.generate_exp_id:
        exp_id = GetTimeString()
        db.Insert('tecan_experiments', [exp_id, -1, "Automatically generated"])
        print "Experiment ID: " + exp_id
    
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
        db.Insert('tecan_experiments', [exp_id, -1, \
                            "Imported from TAR file on " + GetTimeString()])
        WriteToDatabase(MES, db, exp_id)
    
    elif options.xml_filename:
        if not options.plate_id:
            print "Plate ID not supplied, but is mandatory when importing an XML"
            sys.exit(-1)
        if not os.path.exists(options.xml_filename):
            print "File not found: " + options.xml_filename
            sys.exit(-1)
        
        print "Importing from file: " + options.xml_filename
        MES = CollectDataFromSingleFile(options.xml_filename, options.plate_id)

        exp_id = GetCurrentExperimentID(db)
        print "Experiment ID: " + exp_id
        
        WriteToDatabase(MES, db, exp_id)
   
    db.Commit()
    print "Done!"
    
if __name__ == '__main__':
    main()