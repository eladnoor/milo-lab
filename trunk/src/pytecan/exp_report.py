import sys
import numpy as np
from argparse import ArgumentParser
from toolbox.database import MySQLDatabase
import smtplib
 
SMTP_SERVER = 'smtp.gmail.com'
SMTP_PORT = 587

password = 'alphadelta'
 
sender = 'niv.anto@gmail.com'
recipient = 'niv.anto.projects@gmail.com'
subject = 'Test'

 


def MakeOpts():
    """Returns an OptionParser object with all the default options."""
    parser = ArgumentParser()

    parser.add_argument('worklist', nargs=1,
                        help='the path to the worklist that will be written')
    parser.add_argument('-o', '--host', dest='host', default='hldbv02',
                        help='the hostname for the MySQL database')
    parser.add_argument('-d', '--debug', action='store_true', default=False,
                        help='debug mode, store results in dummy DB')
    parser.add_argument("-p", "--num_plates", default=None, type=int, required=True,
                        help="The number of plates in the experiment")
    parser.add_argument("-i", "--iteration", default=None, type=int, required=True,
                        help="The iteration number in the robot script")
    parser.add_argument('-r', '--reading_label', dest='reading_label', default='OD600',
                        help='the label of the measurements used for turbidity tracking')
    return parser

def GetLastPlate(db, plate_num, reading_label):
    max_time = None
    for res in db.Execute('SELECT exp_id, time FROM tecan_readings WHERE reading_label="%s" AND plate=%d ORDER BY time DESC LIMIT 1'
                          % (reading_label, plate_num)):
        exp_id, max_time = res
        break
    
    if max_time is None:
        print "Error in database"
        sys.exit(-1)

    print "Exp ID: %s, Plate: %d, Time: %d" % (exp_id, plate_num, max_time)
    return exp_id, max_time

def GetMeasuredData(db, exp_id, time, plate, reading_label):
    data = np.zeros((8, 12))
    for res in db.Execute('SELECT row, col, measurement FROM tecan_readings WHERE exp_id="%s" AND plate=%d AND time=%d AND reading_label="%s"'
                          % (exp_id, plate, time, reading_label)):
        row, col, measurement = res
        data[row, col] = measurement
    return data
   
def main():

    options = MakeOpts().parse_args()

    db = MySQLDatabase(host=options.host, user='ronm', port=3306,
                       passwd='a1a1a1', db='tecan')
    
    plate_id = options.iteration % options.num_plates

    exp_id, max_time = GetLastPlate(db, plate_id, options.reading_label)
    subject = "%s --> plate : %d " % (exp_id,plate_id)
    data = GetMeasuredData(db, exp_id, max_time, plate_id, options.reading_label)
    
    csv ='' 
    string=''
    for c in range(12):  
        string += "<br>"
        for r in range(8):
            #string += data[r,c].astype('|S6')
            string += "%d , %d --> %s"  % (r+1,c+1,data[r,c].astype('|S6'))
            string += '<br>' 
            print "%d , %d --> %f"  % (r+1,c+1,data[r,c]) 

    for r in range(8):  
        csv += "<br>"
        for c in range(12):
            csv += data[r,c].astype('|S6')
            if (c < 11) :
                csv += ' ,'
            print "%d , %d --> %f"  % (r+1,c+1,data[r,c])              
            
            
            

    #f = open(options.worklist[0], 'w')
    #f.write(csv)
    #f.close() 
    
    print string
    body = string
    body += '<br> <br>'
    body += csv
    headers = ["From: " + sender,
    "Subject: " + subject,
    "To: " + recipient,
    "MIME-Version: 1.0",
    "Content-Type: text/html"]
    
    headers = "\r\n".join(headers)
    session = smtplib.SMTP(SMTP_SERVER, SMTP_PORT)
    session.ehlo()
    session.starttls()
    session.login(sender, password)
    session.sendmail(sender, recipient, headers + "\r\n\r\n" + body)
    session.quit()
    
    
    
    
    print "Done!"
    sys.exit(1)
   
if __name__ == '__main__':
    main()
