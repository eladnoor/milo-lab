import time, calendar, os, csv
from xml.etree.ElementTree import ElementTree
import tarfile
import pylab

fmt = "%Y-%m-%dT%H:%M:%S"

def GetPlateFiles(tar_fname, number_of_plates):
    PL = {}
    L = []
    
    tar = tarfile.open(tar_fname)
    
    for f in tar.getnames():
        if f[-3:] == 'xml':
            ts = f[-23:-3] # take only the postfix of the filename - a timestamp
            ff = tar.extractfile(f)
            L.append((ts, ff))
    
    for i, (_t, f) in enumerate(sorted(L)):    
        PL.setdefault(i % number_of_plates, []).append(f)
    
    return PL    

def ParseReaderFile(fname):
    xml_reader = ElementTree()
    xml_reader.parse(fname)
    well = (0, 0)
    reading_label = ""
    measurement = ""
    TSEC = None
    DATA = {}
    for e in xml_reader.getiterator():
        if e.tag == 'Section':
            reading_label = e.attrib['Name']
            TIME = e.attrib['Time_Start']
            TIME = TIME[:19]
            TS = time.strptime(TIME, fmt)
            if not TSEC:
                TSEC = calendar.timegm(TS)
            DATA[reading_label] = {}
            
        if e.tag == 'Well':
            W = e.attrib['Pos']
            well_row = ord(W[0]) - ord('A')
            well_col = int(W[1:]) - 1
            well = (well_row, well_col)

        if e.tag == 'Multiple':
            if e.attrib['MRW_Position'] == 'Mean':
                measurement = e.text
                DATA[reading_label][well] = float(measurement)

        if e.tag == 'Single':
            measurement = e.text
            if measurement == "OVER":
                DATA[reading_label][well] = pylab.nan
            else:
                DATA[reading_label][well] = float(measurement)

    return DATA, TSEC

def CollectData(tar_fname, number_of_plates):
    PL = GetPlateFiles(tar_fname, number_of_plates)
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

def WriteCSV(MES, csv_fname):
    """
        Write the data into a directory, each reading-label in its own CSV file.
        The columns of the CSV file are: reading-label, plate, well, time, measurement.
        The rows is ordered according to these columns.
    """
    for reading_label in sorted(MES.keys()):
        csv_writer = csv.writer(open(csv_fname, 'well'))
        csv_writer.writerow(['reading label', 'plate', 'row', 'col', 
                             'time', 'measurement'])
        for plate_id in sorted(MES[reading_label].keys()):
            for well in sorted(MES[reading_label][plate_id].keys()):
                for time in sorted(MES[reading_label][plate_id][well].keys()):
                    csv_writer.writerow([reading_label, plate_id, well[0], 
                                         well[1], time, 
                                         MES[reading_label][plate_id][well][time]])

if __name__ == "__main__":
    tar_fname = "../data/tecan/PL6-96.tar.gz"
    csv_fname = "../data/tecan/PL6-96.csv"
    MES = CollectData(tar_fname, number_of_plates=4)
    WriteCSV(MES, csv_fname)
    
    