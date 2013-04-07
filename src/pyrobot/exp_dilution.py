from argparse import ArgumentParser
import tecan
import sys
from database import MySQLDatabase
import numpy as np

def Header():
    """
    Creates a worklist file and writes the first init command (tip size to be used in the fllowing commands)
    """
    header = [] 
    header += ['S;1']
    header +=  ['Vector("Liconic","69","1",0,1,0,1,0,0);']
    header +=  ['Vector("Grid40","40","3",0,1,2,0,0,0);']
    
    return '\n'.join(header)
    
def Footer():
    """
    Creates a worklist file and writes the first init command (tip size to be used in the fllowing commands)
    """
    footer = []
    footer += ['Vector("Grid40","40","3",0,0,0,1,0,0);']
    footer +=  ['Vector("Liconic","69","1",0,1,2,0,0,0);']
     
    return '\n'.join(footer)

    
def Comm(x,labware,row,col,vol,liq):
    """
    Append a command into the worklist file
    Params --> x='A' for aspiarte or x='D' for dispense (char)
               labware = labware position tag on the robot table (string)
               row , col = well cordintates --> converted into 1..96 well position via (col-1)*8+row (int)
               vol = volume in ul  (int)
               liq = liquid class (string)
    """
    pos = (col)*8+row + 1
    return '%s;%s;;;%d;;%d;%s;;;' % (x,labware,pos,vol,liq)
    
def Tip():
    """
    Append a change tup command to worklist file
    """
    return 'W;'
    
def UserPrompt(msg):
    """
    Append a UserPrompt evoke command to the worklist file
    text will be displayed to evoware user while executing this worklist command
    """
    return 'B;UserPrompt("%s",1,-5);' % (msg)

def MakeOpts():
    """Returns an OptionParser object with all the default options."""
    parser = ArgumentParser()

    parser.add_argument("-o", "--host", dest="host", default="hldbv02",
                        help="The hostname for the MySQL database")
    parser.add_argument('-d', '--debug', action='store_true', default=False,
                        help='debug mode, store results in dummy DB')
    parser.add_argument('-w', '--worklist', dest='worklist',
                        help='the path to the worklist that will be written')
    parser.add_argument('-t', '--threshold', dest='threshold', default=0.2, type=float,
                        help='the OD threshold for dilution')
    parser.add_argument('-v', '--volume', dest='vol', default=15, type=int,
                        help='volume for diluation in ul')
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

def get_dilution_rows(db, exp_id, plate, time):
    """
        returns a vector (length 12) of the current rows that should be
        checked for dilution. If they don't exist returns zeros.
    """

    res = db.Execute('SELECT COUNT(*) FROM exp_dilution_columns WHERE exp_id="%s" AND plate=%d'%  (exp_id, plate))
    if res[0][0] == 0:
        for col in xrange(12):
            db.Execute('INSERT INTO exp_dilution_columns(exp_id, plate, col, row, time) VALUES ("%s", %d, %d, 0, %d)' %
                       (exp_id, plate, col, time))

    dilution_rows = [0] * 12
   
    # if the exp_id doesn't exist this loop will not be skipped and the result
    # will be only 0s
    print 'SELECT col, row FROM exp_dilution_columns WHERE exp_id="%s" AND plate=%d' % (exp_id, plate)
    for res in db.Execute('SELECT col, max(row) FROM exp_dilution_columns WHERE exp_id="%s" AND plate=%d GROUP BY col'
                          % (exp_id, plate)):
        col, row = res
        dilution_rows[col] = row
    return dilution_rows

def increment_row(db, exp_id, plate, col, row, time):
    db.Execute('INSERT INTO exp_dilution_columns(exp_id, plate, col, row, time)  VALUES ("%s", %d, %d, %d, %d)' %
               (exp_id, plate, col, row, time))

def main():

    options = MakeOpts().parse_args()
    VOL = options.vol
    LABWARE = 'LABWARE' 
    LIQ = 'LIQUID_CLASS'
    
    # We should also state which directory where the evoware could find the worklist file

    db = MySQLDatabase(host=options.host, user='ronm', port=3306,
                       passwd='a1a1a1', db='tecan')
    
    exp_id, plate, max_time = get_last_plate(db)
    dilution_rows = get_dilution_rows(db, exp_id, plate, max_time)
    
    print "dilution_rows: ", dilution_rows
    
    data = np.zeros((8, 12))
    for res in db.Execute('SELECT row, col, measurement FROM tecan_readings WHERE time=%d AND reading_label="OD600"' % max_time):
        row, col, measurement = res
        data[row, col] = measurement

    worklist = []
    for col, row in enumerate(dilution_rows):
        meas = data[row, col]
        print col, row, meas
        if (meas > options.threshold) and (row < 7):
            msg = "OD = %f --> dilute cell %s%d into cell %s%d" % (meas, chr(ord('A') + row), col+1, chr(ord('A') + row + 1), col+1)
            print msg
            worklist += [UserPrompt(msg)]
            worklist += [Comm('A',LABWARE,row,col,VOL,LIQ)]
            worklist += [Comm('D',LABWARE,row+1,col,VOL,LIQ)]
            #labware,volume and liquid_class would be hard coded for now ...
            worklist += [Tip()]
            increment_row(db, exp_id, plate, col, row+1, max_time)
        
    if len(worklist) > 0:
        worklist = [Header()] + worklist + [Footer()]

    f = open(options.worklist, 'w')
    f.write('\n'.join(worklist))
    f.close()
    db.Commit()
    print "Done!"
    sys.exit(0)
   
if __name__ == '__main__':
    main()
