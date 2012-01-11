#!/usr/bin/python

import logging
import pylab
import numpy
from scipy import stats



def ZeroMatRows(mat):
    Nr, Nc = mat.shape
    mins = numpy.min(mat, 1).reshape(Nr, 1)
    repeat_mins = numpy.repeat(mins, Nc, 1)
    return mat - repeat_mins


class CultureReporterFilterer(object):

    def __init__(self,
                 min_culture_level,
                 min_reporter_level):
        self.min_culture_level = min_culture_level
        self.min_reporter_level = min_reporter_level
    
    def Filter(self, levels_mat, reporter_mat,
               times_mat, labels):
        # Zero the levels
        Nr, Nc = levels_mat.shape
        zeroed_levels = ZeroMatRows(levels_mat)
        zeroed_reporter = ZeroMatRows(reporter_mat)
        
        max_levels = numpy.max(zeroed_levels, 1)
        max_reporter = numpy.max(zeroed_reporter, 1)
        
        high_levels = pylab.find(max_levels >= self.min_culture_level)
        high_reporters = pylab.find(max_reporter >= self.min_reporter_level)
        
        to_keep = set(high_levels)
        to_keep.update(high_reporters)
        to_keep = numpy.array(sorted(to_keep))
        
        return to_keep
        
        
class CultureShifter(object):

    def __init__(self):
        pass
        
    def ShiftCultureTimes(self, levels_mat, times_mat):
        # Zero the levels and times
        Nr, Nc = levels_mat.shape
        zeroed_levels = ZeroMatRows(levels_mat)
        zeroed_times = ZeroMatRows(times_mat)

        # Scale the levels to 1 max        
        maxes = numpy.max(zeroed_levels, 1)
        left_mat = numpy.diag(1/maxes)

        scaled = numpy.dot(left_mat, zeroed_levels)
        min_is = numpy.argmin(numpy.abs(scaled - 0.2), 1)
        
        # Set all cells to cross 20% culture level at time 0.
        shifted_times = zeroed_times.copy()
        for row, col in enumerate(min_is):
            tval = zeroed_times[row, col]
            shifted_times[row, :] -= tval
        
        # Rezero the times.
        min_times = numpy.min(shifted_times)
        shifted_times -= min_times
        
        return shifted_times, scaled
        

class ReporterBackgroundSubtracter(object):
    
    def __init__(self, background_label='BACKGROUND'):
        self.background_label = background_label
        
    def SubtractBackground(self, reporter_mat, levels_mat, labels):
        bg_index = pylab.find(labels == self.background_label)
        assert bg_index.any()
        
        # Assume there's only 1...
        bg_index = bg_index[0]
        bg_reporter = reporter_mat[bg_index, :]
        bg_levels = levels_mat[bg_index, :]
        
        Nr, Nc = reporter_mat.shape
        bg_reporter_mat = numpy.zeros((Nr, Nc))
                
        for i in xrange(Nr):
            for j in xrange(Nc):
                level = levels_mat[i, j]
                bg_min = numpy.abs(level - bg_levels)
                bg_min_idx = numpy.argsort(bg_min)
                closest_bg_i = None
                for k in bg_min_idx:
                    if bg_min[k] > 0:
                        closest_bg_i = k
                        break

                bg_reporter_mat[i, j] = bg_reporter[closest_bg_i]
        
        return reporter_mat - bg_reporter_mat


class ReporterActivityCalculator(object):
    
    def __init__(self,
                 min_culture_level,
                 max_culture_level,
                 min_reporter_level=0,
                 window_size=7):
        assert window_size % 2 == 1
        self.min_culture_level = min_culture_level
        self.max_culture_level = max_culture_level 
        self.min_reporter_level = min_reporter_level
        self.window_size = window_size

    def CalculateAllActivities(self,
                               culture_levels,
                               reporter_levels,
                               times):
        Nr, Nc = culture_levels.shape
        activities_mat = numpy.zeros((Nr, Nc))
        activities_mat[:,:] = numpy.NAN
        for i in xrange(Nr):
            activities = self.CalculateReporterActivity(
                    culture_levels[i,:], reporter_levels[i,:],
                    times[i,:])
            
            
            if activities is not None:
                activities_mat[i, :] = activities
        
        return activities_mat
    
    def _CalculateActivitiesDelta(self, culture_levels, reporter_levels):
        """Calculates the activities from levels."""
        side_size = (self.window_size - 1) / 2
        
        N = len(culture_levels)
        window_deltas = numpy.zeros(N)
        total_ods = numpy.ones(N)
        for i in xrange(side_size, N - side_size):
            window_deltas[i] = numpy.abs(reporter_levels[i+side_size]
                                         - reporter_levels[i-side_size])
            total_ods[i] = numpy.sum(culture_levels[i-side_size:i+side_size+1])
        activities = window_deltas / total_ods
    
        return activities

    def _MovingAverage(self, array_like):
        """Calculates the activities from levels."""
        side_size = (self.window_size - 1) / 2
        
        N = len(array_like)
        averaged = numpy.zeros(N)
        for i in xrange(side_size, N - side_size):
            averaged[i] = numpy.mean(array_like[i-side_size:i+side_size+1])
    
        return averaged
    
    def CalculateAllMaxActivities(self, culture_levels, activities):
        """Calculates the maximal activities for all conditions.
        
        Args:
            culture_levels: matrix of culture levels, row per well.
            activities: matrix of calculated activities, row per well.
        
        Returns:
            A 1d array of maximal activities, 1 entry per well.
        """
        Nr, Nc = activities.shape
        max_activities = numpy.zeros(Nr)
        max_activities[:] = pylab.NAN
        for i in xrange(Nr):
            ac = activities[i,:]
            cu = culture_levels[i,:]
            max_activities[i] = self._CalculateMaxActivity(cu, ac)
        return max_activities        
    
    def _CalculateMaxActivity(self,
                              culture_levels,
                              activities): 
        """Calculates the maximal promoter activity.
        
        Args:
            culture_levels: the culture levels.
            activities: the activities computed.
        
        Returns:
            The maximal (smoothed) reporter level within the
            allowed range of culture levels.
        """        
        abovemin_i = pylab.find(culture_levels >= self.min_culture_level)
        belowmax_i = pylab.find(culture_levels <= self.max_culture_level)
        if not abovemin_i.any():
            return numpy.NAN
        
        lower_i = abovemin_i[0]
        upper_i = belowmax_i[-1]
        if lower_i >= upper_i:
            return numpy.NAN        
        
        window_size = upper_i - lower_i
        
        activities_window = activities[lower_i:lower_i+window_size]
        diff_activities = self._MovingAverage(activities_window)        
        max_activity = numpy.max(diff_activities)
        if max_activity < 0:
            return numpy.NAN        
        
        return max_activity   

    def CalculateReporterActivity(self,
                                  culture_levels,
                                  reporter_levels,
                                  times):
        """Calculates the promoter activities.
        
        Args:
            culture_levels: the culture levels.
            reporter_levels: the reporter levels.
            times: the times for each measurement.
        
        Returns:
            An array of promoter activities for each time point.
        """
        max_culture_level = numpy.max(culture_levels)
        min_culture_level = numpy.min(culture_levels)
        if max_culture_level < self.min_culture_level:
            return None
        if min_culture_level > self.max_culture_level:
            logging.warning('Something is really wrong!')
            return None
        
        max_promoter_level = numpy.max(reporter_levels)
        if max_promoter_level < self.min_reporter_level:
            return None

        activities = self._CalculateActivitiesDelta(
            culture_levels, reporter_levels)

        abovemin_i = pylab.find(culture_levels >= self.min_culture_level)
        belowmax_i = pylab.find(culture_levels <= self.max_culture_level)
        assert abovemin_i.any()
        
        lower_i = abovemin_i[0]
        upper_i = belowmax_i[-1]
        if lower_i >= upper_i:
            logging.warning('Index collision: (%d, %d)', lower_i, upper_i)
            return None
        
        if lower_i < 0 or upper_i >= len(culture_levels):
            return None        
            
        return activities
    
    
        
    