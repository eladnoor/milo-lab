import time, calendar, csv
from xml.etree.ElementTree import ElementTree
import tarfile
import pylab

fmt = "%Y-%m-%dT%H:%M:%S"

def RowCol2String(row, col):
    return "%s%02d" % (chr(ord('A') + row), col+1)

def GetPlateFiles(tar_fname, number_of_plates=None):
    PL = {}
    L = []
    
    tar = tarfile.open(tar_fname, 'r')
    
    for fname in tar.getnames():
        if fname[-3:] == 'xml':
            ts = fname[-23:-3] # take only the postfix of the filename - a timestamp
            f = tar.extractfile(fname)
            L.append((ts, f))
        elif number_of_plates is None and fname[-3:] == 'txt': # this is a 'tag' file that only the number of plates
            f = tar.extractfile(fname)
            try:
                number_of_plates = int(f.read())
            except TypeError:
                raise Exception("the .txt file indicating the number of plates is corrupt")
            f.close()
    
    if number_of_plates is None:
        raise Exception("cannot determine the number of plates")
    
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
    plate_values = {}
    for e in xml_reader.getiterator():
        if e.tag == 'Section':
            reading_label = e.attrib['Name']
            TIME = e.attrib['Time_Start']
            TIME = TIME[:19]
            TS = time.strptime(TIME, fmt)
            time_in_sec = calendar.timegm(TS)
            plate_values[reading_label] = {}
            plate_values[reading_label][time_in_sec] = {}
        elif e.tag == 'Well':
            W = e.attrib['Pos']
            well_row = ord(W[0]) - ord('A')
            well_col = int(W[1:]) - 1
            well = (well_row, well_col)
        elif e.tag == 'Multiple':
            if e.attrib['MRW_Position'] == 'Mean':
                measurement = e.text
                plate_values[reading_label][time_in_sec][well] = float(measurement)
        elif e.tag == 'Single':
            measurement = e.text
            if measurement == "OVER":
                plate_values[reading_label][time_in_sec][well] = pylab.nan
            else:
                plate_values[reading_label][time_in_sec][well] = float(measurement)
    return plate_values

def CollectData(tar_fname, number_of_plates=None):
    PL = GetPlateFiles(tar_fname, number_of_plates)
    MES = {}
    
    for plate_id in PL:
        MES[plate_id] = None
        for f in PL[plate_id]:
            plate_values = ParseReaderFile(f)
            if MES[plate_id] == None:
                MES[plate_id] = plate_values
            else:
                for reading_label, label_values in plate_values.iteritems():
                    for time_in_sec, time_values in label_values.iteritems():
                        MES[plate_id][reading_label][time_in_sec] = time_values
    return MES

def CollectDataFromSingleFile(xml_fname, plate_id):
    return {plate_id: ParseReaderFile(xml_fname)}

def WriteCSV(MES, f):
    """
        Write the data into a directory, each reading-label in its own CSV file.
        The columns of the CSV file are: reading-label, plate, well, time, measurement.
        The rows is ordered according to these columns.
    """
    init_time = GetExpInitTime(MES)
    csv_writer = csv.writer(f)
    csv_writer.writerow(['plate', 'reading label', 'row', 'col', 
                         'time', 'measurement'])
    for plate_id, plate_values in sorted(MES.iteritems()):
        for reading_label, label_values in sorted(plate_values.iteritems()):
            for time_in_sec, time_values in sorted(label_values.iteritems()):
                relative_time_in_hr = (time_in_sec - init_time)/3600.0
                for well, value in sorted(time_values.iteritems()):
                    csv_writer.writerow([plate_id, reading_label, 
                        well[0], well[1], "%.3f" % relative_time_in_hr, value])

def GetExpInitTime(MES):
    all_time_values = []
    for plate_values in sorted(MES.values()):
        for label_values in sorted(plate_values.values()):
            all_time_values += label_values.keys()
    if all_time_values:
        return min(all_time_values)
    else:
        raise ValueError("The experiment has no data, cannot find the init time")
                    
def GetCurrentExperimentID(db):
    for row in db.Execute("SELECT max(exp_id) e FROM tecan_experiments"):
        return row[0]
    raise ValueError("Database Error: no experiments present in tecan_experiments table")
            
def GetExpDate(MES):
    init_time = GetExpInitTime(MES)
    return time.strftime('%Y-%m-%d %H:%M:%S', time.gmtime(init_time))
                    
def WriteToDatabase(MES, db, exp_id):
    for plate_id, plate_values in sorted(MES.iteritems()):
        for reading_label, label_values in sorted(plate_values.iteritems()):
            for time_in_sec, time_values in sorted(label_values.iteritems()):
                for well, value in sorted(time_values.iteritems()):
                    db.Insert('tecan_readings', [exp_id, plate_id, 
                        reading_label, well[0], well[1], time_in_sec, value])
    return exp_id

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

    