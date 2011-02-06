import time, calendar, os, csv
from xml.etree.ElementTree import ElementTree

fmt = "%Y-%m-%dT%H:%M:%S"

def GetPlateFiles(dir, number_of_plates):
    PL = {}
    plate_id = 1
    for i in range(1, number_of_plates+1):
        PL[i] = []
    L = []
    
    for f in os.listdir(dir):
        if f[-3:] == 'xml':
            ts = f[-23:]
            ff = os.path.join(dir, f)
            L.append((ts, ff))
    
    for _t, f in sorted(L):    
        PL[plate_id].append(f)
        plate_id += 1
        if plate_id > number_of_plates:
            plate_id = 1
    
    return PL    

def ParseReaderFile(fname):
    T = ElementTree()
    T.parse(fname)
    W = ""
    M = ""
    V = ""
    TIME = ""
    TSEC = None
    DATA = {}
    for e in T.getiterator():
        if e.tag == 'Section':
            M = e.attrib['Name']
            TIME = e.attrib['Time_Start']
            TIME = TIME[:19]
            TS = time.strptime(TIME   ,fmt)
            if not TSEC:
                TSEC = calendar.timegm(TS)
            DATA[M] = {}
            
        if e.tag == 'Well':
            W = e.attrib['Pos']

        if e.tag == 'Multiple':
            if e.attrib['MRW_Position'] == 'Mean':
                V = e.text
                DATA[M][W] = float(V)

        if e.tag == 'Single':
            V = e.text
            DATA[M][W] = float(V)

    return DATA, TSEC

def CollectData(dir, number_of_plates):
    PL = GetPlateFiles(dir, 4)
    TVEC = {}
    MES = {}
    
    for plate_id in PL:
        TVEC[plate_id] = []
        for f in PL[plate_id]:
            DATA, TSEC = ParseReaderFile(f)
            TVEC[plate_id].append(TSEC)
            for m in DATA:
                MES.setdefault(m, {}).setdefault(plate_id, {})
                for w in DATA[m]:
                    MES[m][plate_id].setdefault(w, {})[TSEC] = DATA[m][w]
    return MES

def WriteCSV(MES, dir):
    """
        Write the data into a directory, each reading-label in its own CSV file.
        The columns of the CSV file are: reading-label, plate, well, time, measurement.
        The rows is ordered according to these columns.
    """
    for m in sorted(MES.keys()):
        csv_writer = csv.writer(open(os.path.join(dir, m + ".csv"), 'w'))
        csv_writer.writerow(['reading label', 'plate', 'well', 'time', 'measurement'])
        for p in sorted(MES[m].keys()):
            for w in sorted(MES[m][p].keys()):
                for t in sorted(MES[m][p][w].keys()):
                    csv_writer.writerow([m, p, w, t, MES[m][p][w][t]])

if __name__ == "__main__":
    dir = "../data/tecan/PL6-96"
    MES = CollectData(dir, number_of_plates=4)
    WriteCSV(MES, dir)
    
    