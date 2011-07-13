from optparse import OptionParser
from pytecan.util import WriteToDatabase, CollectData
import sys, os
from toolbox.database import MySQLDatabase
import tkFileDialog

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
        tar_fname = tkFileDialog.askopenfilename(filetypes=
            [('gzip files', '.gz'), ('tar files', '.tar'), ('all files', '.*')])
        if not tar_fname:
            sys.exit(0)
    else:
        tar_fname = options.tar_filename
    
    if not os.path.exists(tar_fname):
        print "File not found: " + tar_fname
        sys.exit(-1)
    
    print "Importing from file: " + options.tar_filename
    
    #print "Importing into database: " + options.sqlite_db_filename
    if options.exp_id:
        print "Experiment ID: " + options.exp_id
    
    MES = CollectData(options.tar_filename)
    db = MySQLDatabase(host='132.77.80.238', user='ronm', 
                       passwd='a1a1a1', db='tecan')
    exp_id = WriteToDatabase(MES, db, options.exp_id)
    db.Commit()
    print "Done importing experiment: %s" % exp_id
    
    #for row in db.Execute('SELECT exp_id, plate, count(*) from tecan_readings '
    #                      'GROUP BY exp_id, plate'):
    #    print row
    
if __name__ == '__main__':
    main()