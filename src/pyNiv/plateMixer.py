import sys, os, csv
import numpy as np
from argparse import ArgumentParser
from toolbox.database import MySQLDatabase



def Header():
    """
    Creates a worklist file and writes the first init command (tip size to be used in the fllowing commands)
    """
    header = [] 
    header += ['S;1']
    return header
    
def Footer():
    """
    Creates a worklist file and writes the first init command (tip size to be used in the fllowing commands)
    """
    footer = []
    return footer
    
def Comm(x,labware,pos,vol,liq):
    """
    Append a command into the worklist file
    Params --> x='A' for aspiarte or x='D' for dispense (char)
               labware = labware position tag on the robot table (string)
               row , col = well cordintates --> converted into 1..96 well position via (col-1)*8+row (int)
               vol = volume in ul  (int)
               liq = liquid class (string)
    """
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
    return 'B;UserPrompt("%s",0,10);' % (msg)

def MakeOpts():
    """Returns an OptionParser object with all the default options."""
    parser = ArgumentParser()

    parser.add_argument('worklist', nargs=1,
                        help='the path to the worklist that will be written')
    parser.add_argument('-d', '--debug', action='store_true', default=False,
                        help='debug mode, store results in dummy DB')
    parser.add_argument('-t', '--targetOD', dest='targetOD', default=0.1, type=float,
                        help='the OD targetOD for dilution')
    parser.add_argument('-y', '--yfpOD', dest='yfpOD', default=None, required=False, type=float,
                        help='the OD of control YFP strain in the 50 ml tube')                    
    parser.add_argument('-v', '--volume', dest='vol', default=150, type=int,
                        help='final volume')
    parser.add_argument('-l', '--liquid_class', dest='liquid_class', default='TurbidoClass',
                        help='liquid class to be used in pipetation')
    parser.add_argument('-a', '--path', dest='path_csv', default=None, required=True,
                        help='the name of the CVS file which contains the OD data for both plates to be mixed')
    parser.add_argument('-b', '--blank',  default=0.04, type=float, required=False,
                        help='Blank to substact from OD')                    

    return parser

def ReadPathCsv(options):
    if not os.path.exists(options.path_csv):
        error("cannot find the CVS file with the experiment names: " + options.path_csv)
    
    MEDIA_LABWARE = 'MEDIA'
    CHERRY_LABWARE = 'CHERRY'
    YFP_LABWARE = 'YFP'
    DEST_LABWARE = 'DEST'
    LIQ = 'TurbidoClass'
    
    blank = options.blank
    targetOD = options.targetOD
    TH = 0.1 # Threshold for active well
    WELLS = 97

    workArray = [(0,0,0)]*WELLS
    
    for rowNum, line in enumerate(csv.reader(open(options.path_csv, 'r'))):
        for col in xrange(12):
            pos = (col)*8+rowNum + 1
            cherryOD = float(line[col])
            
            if cherryOD < TH :
                print "**** Error in well [%d][%d] --> OD (%f) lower than TH (%f)****\n " % (rowNum, col, cherryOD , TH)
                workArray[pos] = (0,0,0)
                continue
            
            yfpOD = float(line[col+12])
            cherryVol = int((cherryOD/(targetOD/2))**-1 *  float(options.vol))
            print "Cherry : %s in well [%d,%d] %d --> take %d ul" %( cherryOD, rowNum, pos, col , cherryVol)
            yfpVol = int((yfpOD/(targetOD/2))**-1 *  float(options.vol))
            print "YFP : %s in well [%d,%d] %d --> take %d ul" %( yfpOD, rowNum, pos, col+12 , yfpVol)
            mediaVol = options.vol - ( cherryVol + yfpVol)
            print "Media : Media in well [%d,%d] %d --> take %d ul \n\n" %( rowNum, col ,pos, mediaVol)
            workArray[pos] = (mediaVol,yfpVol,cherryVol)

    worklist = []       
    #write header for worklist file
    
    #add to worklist media pipeting commandes 
    msg = 'Pipeting Media !'
    worklist += [UserPrompt(msg)]
    for well in xrange(1,97):
        print  workArray[well][0]
        worklist += [Comm('A',MEDIA_LABWARE,1,workArray[well][0],LIQ)]
        worklist += [Comm('D',DEST_LABWARE,well,workArray[well][0],LIQ)]
    worklist += [Tip()]
        
    #add to worklist YFP control  strain pipeting commandes
    msg = 'Pipeting YFP control strain !'
    worklist += [UserPrompt(msg)]
    for well in xrange(1,97):
        print workArray[well][1]
        worklist += [Comm('A',YFP_LABWARE,well,workArray[well][1],LIQ)]
        worklist += [Comm('D',DEST_LABWARE,well,workArray[well][1],LIQ)]
    worklist += [Tip()]

    #add to worklist mCherry strains pipeting commandes
    msg = 'Pipeting cherry variants strains !'
    worklist += [UserPrompt(msg)]
    for well in xrange(1,97):
        print workArray[well][1]
        worklist += [Comm('A',CHERRY_LABWARE,well,workArray[well][2],LIQ)]
        worklist += [Comm('D',DEST_LABWARE,well,workArray[well][2],LIQ)]
        worklist += [Tip()]
        
    worklist = Header() + worklist + Footer()
    f = open(options.worklist[0], 'w')
    f.write('\n'.join(worklist))
    f.close()
    print "Done!"
   
            
            
    return 0


      

def error(s):
    print "Error"
    sys.exit(-1)

def main():

    options = MakeOpts().parse_args()
    path_dict = ReadPathCsv(options)
    
    # VOL = options.vol
    # MEDIA_VOL = 150-VOL #volune of fresh media in designated well
    
    # LABWARE = 'GRID40SITE3' 
    # EPNSTAND = 'EpnStand'
    
    # LIQ = options.liquid_class
    
    # # We should also state which directory where the evoware could find the worklist file

    # db = MySQLDatabase(host=options.host, user='ronm', port=3306,
                       # passwd='a1a1a1', db='tecan')
    
    # exp_id_dict, plate_id = read_exp_id_csv(options.exp_id_csv)

    # if options.plate not in exp_id_dict:
        # error('The measured plate (%d) does not have an exp_id in the CSV file' % options.plate)

    # exp_id = exp_id_dict[options.plate]

    # max_time = GetLastPlate(db, exp_id, plate_id, options.reading_label)
    # data = GetMeasuredData(db, exp_id, max_time, plate_id, options.reading_label)
    # path_step_dict = GetPathSteps(db, exp_id, plate_id, max_time, path_dict)
    
    # worklist = []
    # for path_label, path_step in path_step_dict.iteritems():
        # row, col = path_dict[path_label][path_step]
        # meas = data[row, col]
        # print path_label, path_step, col, row, meas
        # if (meas > options.targetOD) and (path_step < len(path_dict[path_label])-1):
            # next_row, next_col = path_dict[path_label][path_step+1]
            # msg = "Current plate is : %d ) %s __ OD = %f --> dilute cell %s%d into cell %s%d" % (options.plate, exp_id, meas, chr(ord('A') + row), col+1, chr(ord('A') + next_row), next_col+1)
            # print msg
            # worklist += [UserPrompt(msg)]
            # worklist += [Comm('A',EPNSTAND,0,0,MEDIA_VOL,LIQ)]
            # worklist += [Comm('D',LABWARE,next_row,next_col,MEDIA_VOL,LIQ)]
            # worklist += [Comm('A',LABWARE,row,col,VOL,LIQ)]
            # worklist += [Comm('D',LABWARE,next_row,next_col,VOL,LIQ)]
            # #labware,volume and liquid_class would be hard coded for now ...
            # worklist += [Tip()]
            # IncrementRow(db, exp_id, plate_id, path_label, path_step+1, max_time)
    
    # db.Commit()
    
    # if len(worklist) == 0:
        # sys.exit(0)
    
    # worklist = Header() + worklist + Footer()
    # f = open(options.worklist[0], 'w')
    # f.write('\n'.join(worklist))
    # f.close()
    # print "Done!"
    # sys.exit(1)
   
if __name__ == '__main__':
    main()
