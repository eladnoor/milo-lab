import csv
from toolbox.database import MySQLDatabase
from math import floor

def LoadCSV2DB(filename):
    db = MySQLDatabase(host='132.77.80.238', user='ronm', 
                       passwd='a1a1a1', db='primero')
    counter = 1
    for row in csv.reader(open(filename, 'r')):  
        loc_id = int(row[1])
        name = row[2]
        seq = row[3]
        seq = seq.replace(" ", "")
        
        if name == "":
            continue
        
        box = int(floor(loc_id / 81)) + 1
        residual = loc_id % 81
        row = chr(int(floor(residual / 9)) + 65)
        residual = residual % 9
        col = residual + 1
        
        db.Insert('primers_primer', [counter, name, box, row, col, seq, '', loc_id])
        counter += 1
        db.Commit()

if __name__ == "__main__":
    LoadCSV2DB("primers.csv")