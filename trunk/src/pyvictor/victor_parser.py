"""
    This program requires the xlrd package, which can be found at: http://pypi.python.org/pypi/xlrd
"""
from xlrd import open_workbook
import numpy as np
import matplotlib.pyplot as plt
import types, os, re, time
from time import strptime, strftime
from tkFileDialog import askopenfilename      

class VictorParser():
    
    def __init__(self):
        self.data = []
        self.measurement_names = []
        self.plate = {}
        
        self.reading_label_map = {'Absorbance @ 600 (1.0s) (A)': 'OD600'}
    
    def parse_excel(self, fname):
        if (not os.path.exists(fname)):
            raise Exception("Cannot locate the Excel file: " + fname)
        wd = open_workbook(fname)
        
        protocol_sheet = wd.sheet_by_index(2)
        measurement_cell = protocol_sheet.row_values(98)[0]
        try:
            measured_on = re.findall('\. ([0-9\/\ :]+)$', measurement_cell)[0]
        except IndexError:
            raise Exception("cannot get measurement date in XLS file: " + fname)
        self.measurement_time = strptime(measured_on, '%d/%m/%Y %H:%M:%S')
        
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
            unused_plate, unused_repeat, well, type = row[0:4]
            
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

    def get_time_string(self):
        return strftime('%Y-%m-%d %H:%M:%S', self.measurement_time)
        
    def get_data(self, m_name, row, col):
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
    
    def show_plate(self, m_name=0):
        plt.figure()
        for r in range(8):
            for c in range(12):
                plt.subplot(8, 12, 1+r*12+c)
                (t, v) = self.get_data(m_name, r, c)
                plt.plot(t, v)
        plt.show()

    def get_growth_rate(self, m_name, row, col, window_size=1.5, start_threshold=0.01, plot_figure=False):
        (time, cell_count) = self.get_data(m_name, row, col)
        return VictorParser.fit_growth(time, cell_count, window_size, start_threshold, plot_figure)

    @staticmethod
    def fit_growth(time, cell_count, window_size, start_threshold=0.01, plot_figure=False):
        
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
    
        # get the window-size in samples
        t_mat = np.matrix(time).T
        c_mat = np.matrix(cell_count).T - min(cell_count)
        if (c_mat[-1, 0] == 0):
            c_mat[-1, 0] = min(c_mat[c_mat > 0])
    
        for i in np.arange(N-1, 0, -1):
            if (c_mat[i-1, 0] <= 0):
                c_mat[i-1, 0] = c_mat[i, 0]
    
        c_mat = np.log(c_mat)
        
        res_mat = np.zeros((N,3)) # columns are: slope, offset, error
        for i in range(N):
            try:
                # calculate the indices covered by the window
                i_range = get_frame_range(time, i, window_size)
                x = np.hstack([t_mat[i_range, 0], np.ones((len(i_range), 1))])
                y = c_mat[i_range, 0]
                if (min(np.exp(y)) < start_threshold): # the measurements are still too low to use (because of noise)
                    raise ValueError()
                (a, residues) = np.linalg.lstsq(x, y)[0:2]
                res_mat[i, 0] = a[0]
                res_mat[i, 1] = a[1]
                res_mat[i, 2] = residues
            except ValueError:
                pass
    
        max_i = res_mat[:,0].argmax()
        
        if (plot_figure):
            plt.hold(True)
            plt.plot(time, cell_count-min(cell_count))
            plt.plot(time, res_mat[:,0])
            plt.plot([0, time.max()], [start_threshold, start_threshold], 'r--')
            i_range = get_frame_range(time, max_i, window_size)
            
            x = np.hstack([t_mat[i_range, 0], np.ones((len(i_range), 1))])
            y = x * np.matrix(res_mat[max_i, 0:2]).T
            
            plt.plot(x[:,0], np.exp(y), 'k:', linewidth=4)
            
            #plot(time, errors / errors.max())
            plt.yscale('log')
            #legend(['OD', 'growth rate', 'error'])
            plt.legend(['OD', 'growth rate', 'threshold', 'fit'])
        
        return res_mat[max_i, 0]
    
    @staticmethod
    def fit_growth2(time, cell_count, plot_figure=False):
        def peval(t, p):
            (gr, y_min, y_max, t50) = p
            return y_min + (y_max-y_min)/(1 + np.exp(-gr*(t-t50)))

        def residuals(p, y, t):  
            err = y - peval(t, p) 
            return err
        
        from scipy.optimize import leastsq
        p0 = (1, min(cell_count), max(cell_count), np.mean(time))
        plsq = leastsq(residuals, p0, args=(cell_count, time))[0]
        if (plot_figure):
            plt.plot(time, cell_count, '+')
            plt.plot(time, [peval(t, plsq) for t in time], 'r-')
        
        return plsq[0]

    def write_to_database(self, db, exp_id, plate_id=0):
        for m_name in sorted(self.plate.keys()):
            for (row, col) in sorted(self.plate[m_name].keys()):
                for t, v in self.plate[m_name][(row, col)]:
                    time_in_sec = time.mktime(self.measurement_time) + 3600.0*t
                    db_row = [exp_id, plate_id, m_name, row, col, time_in_sec, v]
                    #print ', '.join(['%s' % x for x in db_row])
                    db.Insert('tecan_readings', db_row)

def callback():
    askopenfilename() 
        
if (__name__ == "__main__"):
    vp = VictorParser()
    fname = askopenfilename(filetypes=[("excel", "*.xls"), ("All files", "*")])
    vp.parse_excel(fname)
    vp.show_plate()
