from argparse import ArgumentParser
import tecan
import sys
from database import MySQLDatabase
import numpy as np

def MakeOpts():
    """Returns an OptionParser object with all the default options."""
    parser = ArgumentParser()

    parser.add_argument("-o", "--host", dest="host", default="hldbv02",
                        help="The hostname for the MySQL database")
    parser.add_argument('-d', '--debug', action='store_true', default=False,
                        help='debug mode, store results in dummy DB')
    parser.add_argument('-w', '--worklist', dest='worklist',
                        help='the path to the worklist that will be written')
    parser.add_argument('-t', '--threshold', dest='threshold', default=0.5,
                        help='the OD threshold for dilution')
    
    return parser

def get_last_plate(db):
    max_time = None
    for res in db.Execute('SELECT exp_id, plate, time FROM tecan_readings WHERE reading_label="OD600" ORDER BY time DESC LIMIT 1'):
        exp_id, plate, max_time = res
        break
    
    if max_time is None:
        print "Error in database"
        sys.exit(-1)

    print "Exp ID: %s, Plate: %d, Time: %d" % (exp_id, plate, max_time)
    return exp_id, plate, max_time

def get_dilution_rows(db, exp_id, plate):
    """
        returns a vector (length 12) of the current rows that should be
        checked for dilution. If they don't exist returns zeros.
    """

    res = db.Execute('SELECT COUNT(*) FROM exp_dilution_columns WHERE exp_id="%s" AND plate=%d'%  (exp_id, plate))
    if res[0][0] == 0:
        for col in xrange(12):
            db.Execute('INSERT INTO exp_dilution_columns(exp_id, plate, col, row) VALUES ("%s", %d, %d, 0)' % (exp_id, plate, col))

    dilution_rows = [0] * 12
   
    # if the exp_id doesn't exist this loop will not be skipped and the result
    # will be only 0s
    print 'SELECT col, row FROM exp_dilution_columns WHERE exp_id="%s" AND plate=%d' % (exp_id, plate)
    for res in db.Execute('SELECT col, row FROM exp_dilution_columns WHERE exp_id="%s" AND plate=%d'
                          % (exp_id, plate)):
        col, row = res
        dilution_rows[col] = row
    return dilution_rows

def write_dilution_rows(db, exp_id, plate, dilution_rows):
    for col, row in enumerate(dilution_rows):
        db.Execute('UPDATE exp_dilution_columns SET row=%d WHERE (exp_id="%s" AND plate=%d AND col=%d)' 
                   % (row, exp_id, plate, col))

def main():
    options = MakeOpts().parse_args()

    db = MySQLDatabase(host=options.host, user='ronm', port=3306,
                       passwd='a1a1a1', db='tecan')
    
    exp_id, plate, max_time = get_last_plate(db)
    dilution_rows = get_dilution_rows(db, exp_id, plate)
    
    print "dilution_rows: ", dilution_rows
    
    data = np.zeros((8, 12))
    for res in db.Execute('SELECT row, col, measurement FROM tecan_readings WHERE time=%d AND reading_label="OD600"' % max_time):
        row, col, measurement = res
        data[row, col] = measurement

    for col, row in enumerate(dilution_rows):
        meas = data[row, col]
        print col, row, meas
        if meas > options.threshold:
            print "dilute cell (%d, %d) into cell (%d, %d)" % (row, col, row+1, col)
            dilution_rows[col] += 1
    
    write_dilution_rows(db, exp_id, plate, dilution_rows)
    
    db.Commit()
    print "Done!"
    sys.exit(0)
   
if __name__ == '__main__':
    main()
