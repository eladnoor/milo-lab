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
    reading_label = None
    time_in_sec = None
    well = (0, 0)
    measurement = None
    DATA = {}
    for e in xml_reader.getiterator():
        if e.tag == 'Section':
            reading_label = e.attrib['Name']
            TIME = e.attrib['Time_Start']
            TIME = TIME[:19]
            TS = time.strptime(TIME, fmt)
            time_in_sec = calendar.timegm(TS)
            DATA[reading_label] = {}
            DATA[reading_label][time_in_sec] = {}
        elif e.tag == 'Well':
            W = e.attrib['Pos']
            well_row = ord(W[0]) - ord('A')
            well_col = int(W[1:]) - 1
            well = (well_row, well_col)
        elif e.tag == 'Multiple':
            if e.attrib['MRW_Position'] == 'Mean':
                measurement = e.text
                DATA[reading_label][time_in_sec][well] = float(measurement)
        elif e.tag == 'Single':
            measurement = e.text
            if measurement == "OVER":
                DATA[reading_label][time_in_sec][well] = pylab.nan
            else:
                DATA[reading_label][time_in_sec][well] = float(measurement)
    return DATA

def CollectData(tar_fname, number_of_plates):
    PL = GetPlateFiles(tar_fname, number_of_plates)
    MES = {}
    
    for plate_id in PL:
        MES[plate_id] = None
        for f in PL[plate_id]:
            DATA = ParseReaderFile(f)
            if MES[plate_id] == None:
                MES[plate_id] = DATA
            else:
                for reading_label, label_DATA in DATA.iteritems():
                    for time_in_sec, time_DATA in label_DATA.iteritems():
                        MES[plate_id][reading_label][time_in_sec] = time_DATA
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

def FitGrowth(time, cell_count, window_size, start_threshold=0.01, plot_figure=False):
    
    def get_frame_range(times, mid_frame, windows_size):
        T = times[mid_frame]
        i_range = []
        for i in range(1, len(times)):
            if (times[i-1] > T - window_size/2.0 and times[i] < T + window_size/2.0):
                i_range.append(i)

        if (len(i_range) < 2): # there are not enough frames to get a good estimation
            raise ValueError()
        return i_range

    N = len(cell_count)
    if (N < window_size):
        raise Exception("The measurement time-series is too short (smaller than the windows-size)")

    t_mat = pylab.matrix(time).T
    
    # normalize the cell_count data by its minimum (
    c_mat = pylab.matrix(cell_count).T - min(cell_count)
    if c_mat[-1, 0] == 0:
        c_mat[-1, 0] = min(c_mat[pylab.find(c_mat > 0)])

    for i in pylab.arange(N-1, 0, -1):
        if c_mat[i-1, 0] <= 0:
            c_mat[i-1, 0] = c_mat[i, 0]

    c_mat = pylab.log(c_mat)
    
    res_mat = pylab.zeros((N, 3)) # columns are: slope, offset, error
    for i in range(N):
        try:
            # calculate the indices covered by the window
            i_range = get_frame_range(time, i, window_size)
            x = pylab.hstack([t_mat[i_range, 0], pylab.ones((len(i_range), 1))])
            y = c_mat[i_range, 0]
            if min(pylab.exp(y)) < start_threshold: # the measurements are still too low to use (because of noise)
                raise ValueError()
            (a, residues) = pylab.lstsq(x, y)[0:2]
            res_mat[i, 0] = a[0]
            res_mat[i, 1] = a[1]
            res_mat[i, 2] = residues
        except ValueError:
            pass

    max_i = res_mat[:,0].argmax()
    
    if plot_figure:
        pylab.hold(True)
        pylab.plot(time, cell_count-min(cell_count))
        pylab.plot(time, res_mat[:,0])
        pylab.plot([0, time.max()], [start_threshold, start_threshold], 'r--')
        i_range = get_frame_range(time, max_i, window_size)
        
        x = pylab.hstack([t_mat[i_range, 0], pylab.ones((len(i_range), 1))])
        y = x * pylab.matrix(res_mat[max_i, 0:2]).T
        
        pylab.plot(x[:,0], pylab.exp(y), 'k:', linewidth=4)
        
        #plot(time, errors / errors.max())
        pylab.yscale('log')
        #legend(['OD', 'growth rate', 'error'])
        pylab.legend(['OD', 'growth rate', 'threshold', 'fit'])
    
    return res_mat[max_i, 0]


if __name__ == "__main__":
    tar_fname = "../data/tecan/PL6-96.tar.gz"
    csv_fname = "../data/tecan/PL6-96.csv"
    MES = CollectData(tar_fname, number_of_plates=4)
    WriteCSV(MES, csv_fname)
    
    