"""
    This program requires the xlrd package, which can be found at: http://pypi.python.org/pypi/xlrd
"""
from xlrd import open_workbook
import numpy as np
import types, re, time
from time import strptime, strftime

class VictorParser():
    # all Victor experiments have only one plate, so their ID is 0
    PLATE_ID = 0

    def __init__(self):
        self.data = []
        self.measurement_names = []
        self.plate = {}
        
        self.reading_label_map = {'Absorbance @ 600 (1.0s) (A)': 'OD600',
                                  'Absorbance' : 'OD600'}
    
    @staticmethod
    def GetInfiniteSheetName(wd):
        for worksheet_name in wd.sheet_names():
            # Check A2 cell for "Device: infinite 200Pro" stamp
            cell_value = wd.sheet_by_name(worksheet_name).cell_value(1, 0)
            if re.search('infinite', cell_value):
                return worksheet_name
        return None
    
    def ParseExcel(self, fp):
        wd = open_workbook(file_contents=fp.read())
        self.measurement_time = None
        self.serial_number = None     
    
        # if file is from GCM -> returns sheet_name with data, else return falase. 
        GCM_sheet_name = VictorParser.GetInfiniteSheetName(wd) 
        if GCM_sheet_name is not None:
            # this should not be hard-coded into the script, but we are lazy:

            GCM_worksheet = wd.sheet_by_name(GCM_sheet_name)
            first_data_row = None
            first_data_col = None
            mode = None
            for r in xrange(GCM_worksheet.nrows):
                row = GCM_worksheet.row_values(r)
                if row[0] is None or type(row[0]) == types.FloatType:
                    continue

                if re.search('Start Time:', row[0]):
                    start_time = row[1]
                    self.measurement_time = strptime(start_time, '%m/%d/%Y %H:%M:%S %p')
                    continue
                                   
                if re.search('System', row[0]):
                    self.serial_number = row[4]
                    continue                       
                    
                if re.search('Mode', row[0]):
                    mode = row[4]
                    continue                    
                    
                if re.search('Cycle Nr.', row[0]):
                    first_data_row = r+1
                    for c in xrange(GCM_worksheet.ncols):
                        if row[c] == 'A1':
                            first_data_col = c
                            break
                    continue

            if self.measurement_time is None:
                raise Exception("cannot get Start Time in XLS file")
                                
            if first_data_row is None or first_data_col is None or mode is None: 
                raise Exception("cannot find mode, first data row, or col in XLS file")
            
            m_name = self.reading_label_map[mode]
            self.plate = {m_name : {}}
            
            for col in xrange(first_data_col, GCM_worksheet.ncols):
                # infer well_name which is stored in this col by parsing col header (i.e 'A1')
                header = GCM_worksheet.cell_value(first_data_row-1, col)
                letter = header[0]
                number = header[1:]
                well_row = ord(letter) - ord('A')
                well_col = int(number) - 1
                self.plate[m_name][(well_row, well_col)] = []
                for row in xrange(first_data_row, GCM_worksheet.nrows):
                    OD_val = GCM_worksheet.cell_value(row, col)
                    if not OD_val:
                        continue
                    time = float(GCM_worksheet.cell_value(row, 1))
                    time_in_hours = time / 3600.0
                    self.plate[m_name][(well_row, well_col)].append((time_in_hours, float(OD_val)))
        
        else: # Assume the XLS is a Victor output file
        
            # Get the date and time of the experiment from the Excel sheet
            # called "Protocol".
            protocol_sheet = wd.sheet_by_index(2)
            
            for r in xrange(protocol_sheet.nrows):
                row = protocol_sheet.row_values(r)[0]
                
                serial_num = re.findall('Instrument serial number: \.* ([0-9]+)', row)
                if len(serial_num) > 0:
                    self.serial_number = serial_num[0]
                    continue
                
                measured_on = re.findall('Measured on \.* ([0-9\s\-\/:]+)$', row)
                if len(measured_on) > 0:
                    self.measurement_time = strptime(measured_on[0], '%d/%m/%Y %H:%M:%S')
                    continue
                
                measured_on = re.findall('Measured on \.* ([0-9\s\-\/:]+ [A|P]M)$', row)
                if len(measured_on) > 0:
                    self.measurement_time = strptime(measured_on[0], '%m/%d/%Y %H:%M:%S %p')
                    continue

            if self.measurement_time is None:
                raise Exception("cannot get measurement date in XLS file")
            
            # Get the values of all the measurements from the "List" Excel sheet
            sheet = wd.sheet_by_index(0)
            titles = sheet.row_values(0) # [Plate, Repeat, Well, Type] + [Time, Measurement] * n
            self.measurement_names = []
            for c in range(5, len(titles), 2):
                m_name = str(titles[c])
                if m_name in self.reading_label_map:
                    m_name = self.reading_label_map[m_name]
                if titles[c] not in self.measurement_names:
                    self.measurement_names.append(m_name)
            
            for m_name in self.measurement_names:
                self.plate[m_name] = {}
                for r in range(8):
                    for c in range(12):
                        self.plate[m_name][(r,c)] = []
                
            for r in range(1, sheet.nrows):
                row = sheet.row_values(r)
                _plate, _repeat, well, _type = row[0:4]
                
                well_row = ord(well[0]) - ord('A')
                well_col = int(well[1:]) - 1
                
                for c in range(4, len(row), 2):
                    try:
                        time = float(row[c]) * 24 # convert days to hours
                        value = float(row[c+1])
                        m_name = str(titles[c+1])
                        if m_name in self.reading_label_map:
                            m_name = self.reading_label_map[m_name]
                        self.plate[m_name][(well_row, well_col)].append((time, value))
                    except ValueError:
                        continue

    def GetData(self, m_name, row, col):
        """
            m_name - the type of measurement. Can be an index (int) according to
                     the order the measurements have been done.
        """
        if type(m_name) == types.IntType:
            m_name = self.measurement_names[m_name]
        data_series = self.plate[m_name][(row, col)]
        times = np.array([t for (t, v) in data_series])
        values = np.array([v for (t, v) in data_series])
        return times, values
    
    def WriteToDB(self, db, exp_id, plate_id=0):
        for m_name in sorted(self.plate.keys()):
            for row, col in sorted(self.plate[m_name].keys()):
                for t, v in self.plate[m_name][(row, col)]:
                    time_in_sec = time.mktime(self.measurement_time) + 3600.0*t
                    db_row = [exp_id, plate_id, m_name, row, col, time_in_sec, v]
                    #print ', '.join(['%s' % x for x in db_row])
                    db.Insert('tecan_readings', db_row)

    @staticmethod
    def GetTimeString(t=None):
        t = t or time.localtime()
        return time.strftime('%Y-%m-%d %H:%M:%S', t)

    @staticmethod
    def ImportFileToDB(fp, db, exp_id=None):
        vp = VictorParser()
        vp.ParseExcel(fp)
        
        exp_id = exp_id or VictorParser.GetTimeString(vp.measurement_time)   
        print "Experiment ID: " + exp_id
    
        # delete any previous data regarding this exp_id
        db.Execute("DELETE FROM tecan_readings WHERE exp_id='%s'" % exp_id)
        q2 = "DELETE FROM tecan_experiments WHERE exp_id='" + exp_id + \
             "' AND serial_number='" + vp.serial_number + "'"
        print q2
        db.Execute(q2)
        db.Execute("DELETE FROM tecan_plates WHERE exp_id='%s'" % exp_id)
        
        desc = "Imported from XLS file on " + VictorParser.GetTimeString()
        db.Insert('tecan_experiments', [exp_id, vp.serial_number, desc])
        db.Insert('tecan_plates', [exp_id, vp.PLATE_ID, "", None, None])
        vp.WriteToDB(db, exp_id)
        db.Commit()
        
        return exp_id
        
if __name__ == "__main__":
    import sys
    v = VictorParser()
    v.ParseExcel(open(sys.argv[1], 'r'))
    print v.GetData('OD600', 0, 0)
