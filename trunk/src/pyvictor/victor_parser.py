"""
    This program requires the xlrd package, which can be found at: http://pypi.python.org/pypi/xlrd
"""
from xlrd import open_workbook
from pylab import *
import types, sys, os

from Tkinter import *
from tkFileDialog import askopenfilename      
from twisted.web import static

class VictorParser():
    
    def __init__(self):
        self.data = []
        self.measurement_names = []
        self.plate = {}
    
    def parse_excel(self, fname):
        if (not os.path.exists(fname)):
            raise Exception("Cannot locate the Excel file: " + fname)
        wd = open_workbook(fname)
        sheet = wd.sheet_by_index(0)
        titles = sheet.row_values(0) # [Plate, Repeat, Well, Type] + [Time, Measurement] * n
        self.measurement_names = []
        for c in range(5, len(titles), 2):
            if (titles[c] not in self.measurement_names):
                self.measurement_names.append(str(titles[c]))
        
        for m_name in self.measurement_names:
            self.plate[m_name] = {}
            for r in range(8):
                for c in range(12):
                    self.plate[m_name][(r,c)] = []
            
        for r in range(1, sheet.nrows):
            row = sheet.row_values(r)
            [plate, repeat, well, type] = row[0:4]
            
            well_row = ord(well[0]) - ord('A')
            well_col = int(well[1:]) - 1
            
            for c in range(4, len(row), 2):
                try:
                    time = float(row[c]) * 24 # convert days to hours
                    value = float(row[c+1])
                    m_name = str(titles[c+1])
                    self.plate[m_name][(well_row, well_col)].append((time, value))
                except ValueError:
                    continue
        
    def get_data(self, m_name, row, col):
        if (type(m_name) == types.IntType):
            m_name = self.measurement_names[m_name]
        times = []
        values = []
        for (time, value) in self.plate[m_name][(row, col)]:
            times.append(time)
            values.append(value)
        return (array(times), array(values))
    
    def show_plate(self, m_name=0):
        figure()
        for r in range(8):
            for c in range(12):
                subplot(8, 12, 1+r*12+c)
                (t, v) = self.get_data(m_name, r, c)
                plot(t, v)
        show()

    def get_growth_rate(self, m_name, row, col, window_size=1.5, start_threshold=0.01, plot_figure=False):
        (time, cell_count) = self.get_data(m_name, row, col)
        return fit_growth(time, cell_count, window_size, start_threshold, plot_figure)

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
        t_mat = matrix(time).T
        c_mat = matrix(cell_count).T - min(cell_count)
        if (c_mat[-1, 0] == 0):
            c_mat[-1, 0] = min(c_mat[find(c_mat > 0)])
    
        for i in arange(N-1, 0, -1):
            if (c_mat[i-1, 0] <= 0):
                c_mat[i-1, 0] = c_mat[i, 0]
    
        c_mat = log(c_mat)
        
        res_mat = zeros((N,3)) # columns are: slope, offset, error
        for i in range(N):
            try:
                # calculate the indices covered by the window
                i_range = get_frame_range(time, i, window_size)
                x = hstack([t_mat[i_range, 0], ones((len(i_range), 1))])
                y = c_mat[i_range, 0]
                if (min(exp(y)) < start_threshold): # the measurements are still too low to use (because of noise)
                    raise ValueError()
                (a, residues) = lstsq(x, y)[0:2]
                res_mat[i, 0] = a[0]
                res_mat[i, 1] = a[1]
                res_mat[i, 2] = residues
            except ValueError:
                pass
    
        max_i = res_mat[:,0].argmax()
        
        if (plot_figure):
            hold(True)
            plot(time, cell_count-min(cell_count))
            plot(time, res_mat[:,0])
            plot([0, time.max()], [start_threshold, start_threshold], 'r--')
            i_range = get_frame_range(time, max_i, window_size)
            
            x = hstack([t_mat[i_range, 0], ones((len(i_range), 1))])
            y = x * matrix(res_mat[max_i, 0:2]).T
            
            plot(x[:,0], exp(y), 'k:', linewidth=4)
            
            #plot(time, errors / errors.max())
            yscale('log')
            #legend(['OD', 'growth rate', 'error'])
            legend(['OD', 'growth rate', 'threshold', 'fit'])
        
        return res_mat[max_i, 0]
    
    @staticmethod
    def fit_growth2(time, cell_count, plot_figure=False):
        def peval(t, p):
            (gr, y_min, y_max, t50) = p
            return y_min + (y_max-y_min)/(1 + exp(-gr*(t-t50)))

        def residuals(p, y, t):  
            err = y - peval(t, p) 
            return err
        
        from scipy.optimize import leastsq
        p0 = (1, min(cell_count), max(cell_count), mean(time))
        plsq = leastsq(residuals, p0, args=(cell_count, time))[0]
        if (plot_figure):
            plot(time, cell_count, '+')
            plot(time, [peval(t, plsq) for t in time], 'r-')
        
        return plsq[0]

def callback():
    askopenfilename() 
        
if (__name__ == "__main__"):
    vp = VictorParser()
    fname = askopenfilename(filetypes=[("excel", "*.xls"), ("All files", "*")])
    vp.parse_excel(fname)
    vp.show_plate()
