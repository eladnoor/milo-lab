#!/usr/bin/python

import pylab
import numpy
import scipy.optimize


class GrowthCalculator(object):
    
    def __init__(self):
        pass
    
    @staticmethod
    def NormalizeTimes(times):
        return times - min(times)

    def CalculateGrowthInternal(self, times, levels):
        """Internal implementation of CalculateGrowth.
        
        To be implemented by subclasses.
        
        Returns:
            A two-tuple (growth rate, stationary level)
        """
        raise NotImplementedError
    
    def CalculateGrowth(self, times, levels):
        """Calculate the growth rate and stationary level.
        
        Args:
            times: a 1d numpy.array of timestamps in seconds.
            levels: a 1d numpy.array of levels for each timestamp.
            baseline_levels: a 1d numpy.array of levels per timestamp
                in a baseline condition.
        """
        new_times = self.NormalizeTimes(times)        
        return self.CalculateGrowthInternal(new_times, levels)
    
    def CalculatePlateGrowth(self, plate, reading_label):
        """Calculate the growth rates for a whole plate."""
        times, readings, labels = plate.SelectReading(reading_label)
        
        n = labels.size
        rates = {}
        stationaries = {}
        for i in xrange(n):
            label = labels[i]
            
            rate, stationary = self.CalculateGrowth(times[i,:], readings[i,:])
            scaled_rate = rate * 60 * 60
            
            rates.setdefault(label, []).append(scaled_rate)
            stationaries.setdefault(label, []).append(stationary)
        
        return rates, stationaries

    
    def CalculateStationary(self, times, levels):
        """Return the stationary level."""
        return self.CalculateGrowth(times, levels)[1]

    def CalculateGrowthRate(self, times, levels):
        """Return the stationary level."""
        return self.CalculateGrowth(times, levels)[0]
    
    
class SlidingWindowGrowthCalculator(GrowthCalculator):
    
    def __init__(self,
                 window_size=5,
                 minimum_level=0.01,
                 maximum_level=0.3):
        GrowthCalculator.__init__(self)
        self.window_size = window_size
        self.minimum_level = minimum_level
        self.maximum_level = maximum_level

    def CalculateRates(self, times, levels):
        N = len(levels)    
        t_mat = pylab.matrix(times).T
        
        # normalize the cell_count data by its minimum
        count_matrix = pylab.matrix(levels).T
        norm_counts = count_matrix - min(levels)
        c_mat = pylab.matrix(norm_counts)
        if c_mat[-1, 0] == 0:
            c_mat[-1, 0] = min(c_mat[pylab.find(c_mat > 0)])
    
        for i in pylab.arange(N-1, 0, -1):
            if c_mat[i-1, 0] <= 0:
                c_mat[i-1, 0] = c_mat[i, 0]
    
        c_mat = pylab.log(c_mat)
        
        res_mat = pylab.zeros((N, 5)) # columns are: slope, offset, error, avg_value, max_value
        for i in xrange(N-self.window_size):
            i_range = range(i, i+self.window_size)
            x = pylab.hstack([t_mat[i_range, 0], pylab.ones((len(i_range), 1))])
            y = c_mat[i_range, 0]
            
            # Measurements in window must all be above the min.
            if min(pylab.exp(y)) < self.minimum_level:
                continue
                        
            (a, residues) = pylab.lstsq(x, y)[0:2]
            res_mat[i, 0] = a[0]
            res_mat[i, 1] = a[1]
            res_mat[i, 2] = residues
            res_mat[i, 3] = pylab.mean(count_matrix[i_range,0])
            res_mat[i, 4] = max(pylab.exp(y))
        
        return res_mat
        
    def CalculateGrowthInternal(self, times, levels):
        res_mat = self.CalculateRates(times, levels)
        t_mat = pylab.matrix(times).T
        count_matrix = pylab.matrix(levels).T
        norm_counts = count_matrix - min(levels)
        allowed_vals = pylab.find(res_mat[:,4] < self.maximum_level)
        max_i = res_mat[allowed_vals,0].argmax()

        abs_res_mat = pylab.array(res_mat)
        abs_res_mat[:,0] = pylab.absolute(res_mat[:,0])
        order = abs_res_mat[:,0].argsort(axis=0)
        stationary_indices = filter(lambda x: x >= max_i, order)
        stationary_indices = pylab.array(filter(lambda x: res_mat[x,3] > 0,
                                                stationary_indices))
        
        stationary_level = 0.0
        if stationary_indices.any():
            stationary_level = res_mat[stationary_indices[0], 3]
        
        """
        pylab.hold(True)
        pylab.plot(times, norm_counts)
        pylab.plot(times, res_mat[:,0])
        pylab.plot([0, times.max()], [self.minimum_level, self.minimum_level], 'r--')
        pylab.plot([0, times.max()], [self.maximum_level, self.maximum_level], 'r--')
        i_range = range(max_i, max_i+self.window_size)
        
        x = pylab.hstack([t_mat[i_range, 0], pylab.ones((len(i_range), 1))])
        y = x * pylab.matrix(res_mat[max_i, 0:2]).T
        pylab.plot(x[:,0], pylab.exp(y), 'k:', linewidth=4)
                
        #pylab.plot([0, max(times)], [stationary_level, stationary_level], 'k-')
        
        pylab.yscale('log')
        pylab.legend(['OD', 'growth rate', 'threshold', 'fit'])
        #, 'stationary'])
        """
        
        return res_mat[max_i, 0], stationary_level


class FitSigmoidGrowthCalculator(GrowthCalculator):
    
    def __init__(self):
        GrowthCalculator.__init__(self)
    
    @staticmethod
    def Sigmoid(params, x):
        x0,y0,c,k = params
        y = c / (1 + numpy.exp(-k*(x-x0))) + y0
        return y
    
    @staticmethod
    def Residuals(params, x, y):
        return y - FitSigmoidGrowthCalculator.Sigmoid(params, x)

    @staticmethod
    def FitSigmoid(x, y):
        p_guess=(numpy.median(x), numpy.median(y),
                 max(y), min(y))
        p, cov, infodict, mesg, ier = scipy.optimize.leastsq(
            FitSigmoidGrowthCalculator.Residuals,
            p_guess, args=(x,y),
            full_output=1)                        
        return p
    
    @staticmethod
    def Midpoint(params):
        x0,y0,c,k = params
        ymid = y0 + (c / 2.0)
        xmid = (-1/k)*numpy.log((c/(ymid-y0)) - 1) + x0
        return xmid, ymid 

    @staticmethod
    def GrowthRate(params):
        x0,y0,c,k = params
        xmid, ymid = FitSigmoidGrowthCalculator.Midpoint(params) 
        e = numpy.exp(-k*(xmid-x0))
        rate = k*c*e / (1+e)**2
        offset = ymid - rate*xmid
        return rate, offset
    
    @staticmethod
    def Stationary(params):
        x0,y0,c,k = params
        return c + y0
    
    @staticmethod
    def SigmoidStr(params):
        x0,y0,c,k = params    
        return '%.2g + %.2g / (1 + exp(%.2g (x - %.2g)))' % (y0, c, -k, x0)
    
    def CalculateGrowthInternal(self, times, levels):
        new_levels = levels / max(levels)
        new_times = times / (100.0*max(times))
        log_levels = numpy.log(new_levels)
        params = self.FitSigmoid(new_times, log_levels)
        
        print 'sigmoid'
        print self.SigmoidStr(params)
                
        pylab.plot(new_times, log_levels, 'b.')
        pylab.plot(new_times, self.Sigmoid(params, new_times), 'g-')
        
        slope, offset = self.GrowthRate(params)
        pylab.plot(new_times, slope*new_times + offset, 'r-')

        stationary = self.Stationary(params)
        pylab.plot(new_times, pylab.ones(len(new_times))*stationary)
        
        return None, None


if __name__ == '__main__':
    
    noisy_vals =[0.0374,0.0363,0.0396,0.0396,0.0394,0.0393,0.0394,0.0394,0.0395,0.0395,
                 0.0397,0.0397,0.0398,0.0395,0.0399,0.0399,0.04,0.0402,0.04,0.0401,0.0402,
                 0.0405,0.0403,0.0406,0.0402,0.0407,0.0409,0.041,0.0411,0.0412,0.0414,0.0418,
                 0.0417,0.0421,0.0422,0.043,0.0429,0.0438,0.044,0.0445,0.0446,0.0453,0.0461,0.0477,
                 0.0485,0.0494,0.0504,0.0519,0.0537,0.0564,0.0577,0.061,0.0645,0.0684,
                 0.0732,0.0791,0.0851,0.0915,0.1012,0.1102,0.1212,0.1333,0.1455,0.1534,
                 0.1615,0.1762,0.1848,0.1962,0.2133,0.2276,0.2386,0.2501,0.2645,0.273,
                 0.2845,0.2947,0.3062,0.3115,0.3189,0.3274,0.3355,0.3401,0.3458,0.3464,
                 0.3494,0.358,0.3508,0.3537,0.3406,0.334,0.3411,0.3301,0.326]
    
    #times = [1.315494075e9,1.31549589e9,1.315497673e9,1.315499473e9,1.315501273e9,1.315503073e9,1.315504873e9,1.315506673e9,1.315508473e9,1.315510273e9,1.315512073e9,1.315513874e9,1.315515672e9,1.315517472e9,1.315519273e9,1.315521073e9,1.315522873e9,1.315524672e9,1.315526473e9,1.315528272e9,1.315530098e9,1.315531873e9,1.315533673e9,1.315535474e9,1.315537274e9,1.315539075e9,1.315540875e9,1.315542673e9,1.315544475e9,1.315546275e9,1.315548073e9,1.315549873e9,1.315551674e9,1.315553493e9,1.315555275e9,1.315557074e9,1.315558875e9,1.315560675e9,1.315562475e9,1.315564276e9,1.315566075e9,1.315567875e9,1.315569675e9,1.315571476e9,1.315573276e9,1.315575075e9,1.315576876e9,1.315578676e9]
    #noisy_vals = [5.0e-3,5.0e-3,5.0e-3,5.0e-3,5.0e-3,5.0e-3,5.0e-3,5.0e-3,5.0e-3,5.0e-3,5.0e-3,5.338541666666655e-3,7.2385416666666536e-3,9.638541666666653e-3,1.1738541666666658e-2,1.7638541666666653e-2,2.2438541666666652e-2,2.8338541666666654e-2,3.013854166666665e-2,3.403854166666665e-2,4.033854166666665e-2,4.693854166666666e-2,5.703854166666666e-2,6.783854166666667e-2,7.773854166666666e-2,9.023854166666664e-2,0.10373854166666666,0.11543854166666664,0.12933854166666664,0.14513854166666665,0.15813854166666666,0.17763854166666665,0.19113854166666666,0.19433854166666664,0.19213854166666666,0.19063854166666666,0.18833854166666664,0.18663854166666666,0.18543854166666665,0.18783854166666664,0.18783854166666664,0.18693854166666665,0.18263854166666665,0.17973854166666664,0.18753854166666664,0.18343854166666665,0.18523854166666665,0.18553854166666664]
    #times = pylab.array(times)
    levels = pylab.array(noisy_vals)
    times = pylab.array(map(float, range(len(noisy_vals))))

    calculator = SlidingWindowGrowthCalculator()
    rate, stationary = calculator.CalculateGrowth(times, levels)
    print rate, stationary
    
    #pylab.plot(times, levels, 'g.')
    #pylab.plot(times[:10], rate*times[:10], 'b-')
    pylab.show()
    
    
        
    