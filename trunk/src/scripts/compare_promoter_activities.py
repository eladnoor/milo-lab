#!/usr/bin/python

"""Calculates promoter activities from a plate."""

import pylab
import sys
import numpy
import scipy.stats as st

from toolbox.database import MySQLDatabase
from toolbox import growth
from toolbox import promoter_activity
from toolbox import util
from toolbox.stats import MeanWithConfidenceInterval
from toolbox.color import ColorMap
from toolbox.plate import Plate96
from scripts import templates

from optparse import OptionParser
from os import path
from matplotlib.font_manager import FontProperties


LEGEND_FONT = FontProperties(size=8)


def MakeOpts():
    """Returns an OptionParser object with all the default options."""
    opt_parser = OptionParser()
    opt_parser.add_option("-e", "--experiment_id",
                          dest="experiment_id",
                          help="experiment ID")
    opt_parser.add_option("-f", "--first_plate_ids",
                          dest="first_plate_ids",
                          help="Plates with the first condition.")
    opt_parser.add_option("-s", "--second_plate_ids",
                          dest="second_plate_ids",
                          help="Plates with the second condition.")
    opt_parser.add_option("-c", "--culture_label",
                          dest="culture_label",
                          default="OD600",
                          help="Culture size label")
    opt_parser.add_option("-r", "--reporter_label",
                          dest="reporter_label",
                          default="GFP",
                          help="Promoter expression label")
    opt_parser.add_option("-b", "--background_label",
                          dest="background_label",
                          help="Label of the background strain/well.")
    opt_parser.add_option("-w", "--window_size",
                          dest="window_size",
                          type='int',
                          default=7,
                          help="Window size for computing the growth rate.")
    opt_parser.add_option("-m", "--min_culture_level",
                          dest="min_culture_level",
                          type='float',
                          default=0.2,
                          help="Minimum culture level to consider valid.")
    opt_parser.add_option("-n", "--min_reporter_level",
                          dest="min_reporter_level",
                          type='float',
                          default=300,
                          help="Minimum reporter level to consider valid.")
    opt_parser.add_option("-l", "--lower_culture_bound",
                          dest="lower_culture_bound",
                          type='float',
                          default=0.05,
                          help="Lower bound on culture window (post processed).")
    opt_parser.add_option("-u", "--upper_culture_bound",
                          dest="upper_culture_bound",
                          type='float',
                          default=0.2,
                          help="Upper bound on culture window (post processed).")
    opt_parser.add_option("-p", "--lower_reporter_bound",
                          dest="lower_reporter_bound",
                          type='float',
                          default=10,
                          help="Upper bound on reporter (post processed).")
    opt_parser.add_option("-o", "--output_dir",
                          dest="output_dir",
                          default='../res/activity/',
                          help="Directory to write output to.")
    opt_parser.add_option("-g", "--labels_to_ignore",
                          dest="labels_to_ignore",
                          default="BLANK",
                          help="Labels to ignore, comma separated.")
    return opt_parser


def GetMeansAndErrors(activities, labels):
    mean_and_err_rates = [MeanWithConfidenceInterval(activities[l])
                          for l in labels]
    means = [t[0] for t in mean_and_err_rates]
    errs = [t[1] for t in mean_and_err_rates]
    return means, errs

def UpdateActivitiesDict(activities, new_activities):
    for k, val in new_activities.iteritems():
        activities.setdefault(k, []).extend(val)
    return activities

def GetGeneTextColor(label):
    if 'ED' in label:
        return 'g'
    if 'UPPER GLYC' in label:
        return 'r'
    if 'LOWER GLYC' in label:
        return 'b'
    if 'RIBOSOME' in label:
        return 'c'
    return 'k'


class PlateActivityRunner(object):
    
    def __init__(self,
                 culture_label,
                 reporter_label,
                 filterer,
                 culture_shifter,
                 background_subtracter,
                 activity_calc):
        """Initialize the PlateActivityRunner.
        
        Args:
            plate: the plate.
            culture_label: the culture measurement label.
            reporter_label: the reporter measurement label.
            filterer: a CultureReporterFilterer that filters wells.
            culture_shifter: a CultureShifter object shifts time measurements.
            background_subtracter: ReporterBackgroundSubtracter that 
                subtracts background measurements.
            activity_calc: a ReporterActivityCalculator that calculates
                the activity level.
        """
        self.plates = []
        self.culture_label = culture_label
        self.reporter_label = reporter_label
        self.filterer = filterer
        self.culture_shifter = culture_shifter
        self.reporter_bg_subtracter = background_subtracter
        self.activity_calc = activity_calc
        
    def AddPlate(self, plate):
        self.plates.append(plate)
    
    def GetReadings(self, reading_label):
        times = []
        levels = []
        labels = []
        for plate in self.plates:
            t, lev, lab = plate.SelectReading(reading_label)
            times.append(t)
            levels.append(lev)
            labels.append(lab)
        
        return numpy.vstack(times), numpy.vstack(levels), numpy.hstack(labels)
    
    def Run(self):
        """Runs the analysis. Stores all output in self."""
        # Fetch data from the plate
        culture_times, culture_levels, culture_labels = self.GetReadings(
                self.culture_label)
        unused_times, reporter_levels, unused_labels = self.GetReadings(
                self.reporter_label)
        
        assert reporter_levels.shape[0] == culture_levels.shape[0]
        
        self.raw_times = culture_times
        self.raw_labels = culture_labels
        self.raw_culture_levels = culture_levels
        self.raw_reporter_levels = reporter_levels 
        
        # Filter data from the plate and record what we filtered.
        self.indices_kept = self.filterer.Filter(
            self.raw_culture_levels, self.raw_reporter_levels,
            self.raw_times, self.raw_labels)
        self.filtered_labels = self.raw_labels[self.indices_kept]
        kept_set = set(self.indices_kept)
        self.indices_removed = numpy.array([i for i in xrange(len(self.raw_labels))
                                            if i not in kept_set])
        self.labels_removed = numpy.array([])
        if self.indices_removed.any():
            self.labels_removed = self.raw_labels[self.indices_removed]
        
        self.filtered_times = self.raw_times[self.indices_kept, :]
        self.filtered_culture_levels = self.raw_culture_levels[self.indices_kept, :]
        self.filtered_reporter_levels = self.raw_reporter_levels[self.indices_kept, :]
        
        # Shift the culture times to align the growth in time.
        self.shifted_times, self.scaled_culture_levels = self.culture_shifter.ShiftCultureTimes(
            self.filtered_culture_levels, self.filtered_times)
        
        # Subtract the background from the reporter levels.
        self.filtered_reporter_levels_no_bg = self.reporter_bg_subtracter.SubtractBackground(
            self.filtered_reporter_levels, self.filtered_culture_levels,
            self.filtered_labels)
        
        # Compute the activities of the strains over time.
        self.filtered_activities = self.activity_calc.CalculateAllActivities(
            self.filtered_culture_levels, self.filtered_reporter_levels_no_bg,
            self.filtered_times)
        
        # Compute the smoothed activities of strains over time.
        self.smooth_filtered_activities = self.activity_calc.SmoothAllActivities(
            self.filtered_activities)
        
        # Compute the maximal activities of each strain.
        self.filtered_max_activities = self.activity_calc.CalculateAllMaxActivities(
            self.filtered_culture_levels, self.filtered_activities)
        
        # Make a dict of the maximal activities.
        self.filtered_max_activities_dict = {}
        for k, v in zip(self.filtered_labels, self.filtered_max_activities):
            self.filtered_max_activities_dict.setdefault(k, []).append(v)


class StrainData(object):
    
    def __init__(self, label):
        self.label = label
        self.conditions_to_raw_culture_levels = {}
        self.conditions_to_raw_reporter_levels = {}
        self.conditions_to_raw_times = {}
        self.conditions_to_culture_levels = {}
        self.conditions_to_reporter_levels = {}
        self.conditions_to_times = {}
        self.conditions_to_activities = {}
        self.conditions_to_smooth_activities = {}
        
        self.slug_label = util.slugify(self.label)
        self.raw_levels_fname = '%s_raw_levels.png' % self.slug_label
        self.levels_fname = '%s_levels.png' % self.slug_label
        self.vs_bg_fname = '%s_vs_bg.png' % self.slug_label
        self.activity_fname = '%s_activity.png' % self.slug_label
    
    def AddConditionRawCultureLevels(self, condition, culture_levels):
        self.conditions_to_raw_culture_levels.setdefault(condition, []).append(culture_levels)
    
    def AddConditionRawReporterLevels(self, condition, reporter_levels):
        self.conditions_to_raw_reporter_levels.setdefault(condition, []).append(reporter_levels)
    
    def AddConditionRawTimes(self, condition, times):
        self.conditions_to_raw_times.setdefault(condition, []).append(times)
    
    def AddConditionCultureLevels(self, condition, culture_levels):
        self.conditions_to_culture_levels.setdefault(condition, []).append(culture_levels)
        
    def AddConditionReporterLevels(self, condition, reporter_levels):
        self.conditions_to_reporter_levels.setdefault(condition, []).append(reporter_levels)
    
    def AddConditionTimes(self, condition, times):
        self.conditions_to_times.setdefault(condition, []).append(times)
    
    def AddConditionActivities(self, condition, activities, max_activity):
        t = (activities, max_activity)
        self.conditions_to_activities.setdefault(condition, []).append(t)

    def AddConditionSmoothActivities(self, condition, activities, max_activity):
        t = (activities, max_activity)
        self.conditions_to_smooth_activities.setdefault(condition, []).append(t)

    def GetMeanMaxActivity(self, condition):
        activity_vals = self.conditions_to_activities.get(condition, None)
        if activity_vals is None:
            return numpy.NAN
        
        maxes = [v[1] for v in activity_vals]
        return MeanWithConfidenceInterval(maxes)
    
    def MakeFigures(self, dirname, background_strain):
        colormap = ColorMap(self.conditions_to_raw_times.keys())

        # Raw data figure.
        pylab.figure()
        # Subplot for culture levels  
        pylab.subplot('211')
        pylab.title('Raw Levels for %s' % self.label)
        pylab.xlabel('Time (s)')
        pylab.ylabel('Raw Culture Level')
        for condition, times_list in self.conditions_to_raw_times.iteritems():
            levels_list = self.conditions_to_raw_culture_levels[condition]
            color = colormap[condition]
            for i, (times, levels) in enumerate(zip(times_list, levels_list)):
                pylab.plot(times, levels, color=color, linestyle='-')
        
        # Subplot for reporter levels
        pylab.subplot('212')
        pylab.xlabel('Time (s)')
        pylab.ylabel('Raw Reporter Level')
        for condition, times_list in self.conditions_to_raw_times.iteritems():
            levels_list = self.conditions_to_raw_reporter_levels[condition]
            color = colormap[condition]
            for i, (times, levels) in enumerate(zip(times_list, levels_list)):
                label = None
                if i == 0:
                    label = condition
                    
                pylab.plot(times, levels, color=color, label=label,
                           linestyle='-')
        pylab.legend(loc='upper left', prop=LEGEND_FONT)
        
        # Save the raw data figure.
        pylab.savefig(path.join(dirname, self.raw_levels_fname),
                      format='png')
        
        # Processed data figure.
        pylab.figure()
        # Subplot for culture levels  
        pylab.subplot('211')
        pylab.title('Levels for %s' % self.label)
        pylab.xlabel('Time (s)')
        pylab.ylabel('Culture Level')
        for condition, times_list in self.conditions_to_times.iteritems():
            levels_list = self.conditions_to_culture_levels[condition]
            color = colormap[condition]
            for i, (times, levels) in enumerate(zip(times_list, levels_list)):
                pylab.plot(times, levels, color=color, linestyle='-')
        
        # Subplot for reporter levels
        pylab.subplot('212')
        pylab.xlabel('Time (s)')
        pylab.ylabel('Reporter Level')
        for condition, times_list in self.conditions_to_times.iteritems():
            levels_list = self.conditions_to_reporter_levels[condition]
            color = colormap[condition]
            for i, (times, levels) in enumerate(zip(times_list, levels_list)):
                label = None
                if i == 0:
                    label = condition
                    
                pylab.plot(
                    times, levels, color=color, label=label, linestyle='-')
        pylab.legend(loc='upper left', prop=LEGEND_FONT)
        
        # Save the levels figure.
        pylab.savefig(path.join(dirname, self.levels_fname),
                      format='png')
        
        pylab.figure()
        pylab.subplot('211')
        pylab.title('Activity for %s' % self.label)
        pylab.xlabel('Time (s)')
        pylab.ylabel('Reporter Activity')
        for condition, times_list in self.conditions_to_times.iteritems():
            activities_list = self.conditions_to_activities[condition]
            color = colormap[condition]
            for i, (times, activities_data) in enumerate(zip(times_list, activities_list)):
                activities, max_activity = activities_data
                pylab.plot(times, activities, color=color,
                           linestyle='-')
                pylab.plot(times, [max_activity]*len(times),
                           color=color, linestyle='--')
        
        pylab.subplot('212')
        pylab.xlabel('Time (s)')
        pylab.ylabel('Smoothed Reporter Activity')
        for condition, times_list in self.conditions_to_times.iteritems():
            activities_list = self.conditions_to_smooth_activities[condition]
            color = colormap[condition]
            for i, (times, activities_data) in enumerate(zip(times_list, activities_list)):
                label = None
                max_label = None
                if i == 0:
                    label = condition
                    max_label = 'Max activity in %s' % condition
                
                activities, max_activity = activities_data
                pylab.plot(times, activities, color=color, label=label,
                           linestyle='-')
                pylab.plot(times, [max_activity]*len(times),
                           color=color, label=max_label,
                           linestyle='--')
        
        pylab.legend(loc='lower right', prop=LEGEND_FONT)

        # Save the activity figure.
        pylab.savefig(path.join(dirname, self.activity_fname),
                      format='png')
        
        # Plot reporter levels against their background.
        pylab.figure()
        pylab.title('Raw Reporter vs. Raw Culture (%s and Background)' % self.label)
        pylab.xlabel('Culture Level')
        pylab.ylabel('Reporter Level')

        for condition, reporter_levels_list in self.conditions_to_raw_reporter_levels.iteritems():
            color = colormap[condition]
            culture_levels_list = self.conditions_to_raw_culture_levels[condition]
            bg_reporter_levels_list = background_strain.conditions_to_raw_reporter_levels[condition]
            bg_culture_levels_list = background_strain.conditions_to_raw_culture_levels[condition]
            
            for i, (clevels, rlevels) in enumerate(zip(culture_levels_list, reporter_levels_list)):
                label = None
                if i == 0:
                    label = condition
                pylab.plot(clevels, rlevels, color=color, label=label,
                           linestyle='-')
                
            for i, (clevels, rlevels) in enumerate(zip(bg_culture_levels_list, bg_reporter_levels_list)):
                label = None
                max_label = None
                if i == 0:
                    label = '%s background' % condition
                pylab.plot(clevels, rlevels, color=color, label=label,
                           linestyle='--')
        pylab.legend(loc='upper left', prop=LEGEND_FONT)
        
        # Save the vs background figure.
        pylab.savefig(path.join(dirname, self.vs_bg_fname),
                      format='png')
        
        
class StrainConditionsData(object):
    
    def __init__(self, background_label):
        self.strains = {}
        self.conditions = set()
        self.background_label = background_label
        self.plates = {}
    
    def GetStrainLabels(self):
        return sorted(self.strains.keys())

    def AllStrains(self):
        sorted_keys = self.GetStrainLabels()
        for k in sorted_keys:
            yield self.strains[k]
    all_strains = property(AllStrains)
        
    def AddPlateData(self, condition, plate_runner,
                     ignore_labels=None):
        my_ignore_labels = ignore_labels or set()
        self.conditions.add(condition)
        self.plates.setdefault(condition, []).append(plate_runner)
        
        # Add raw data
        labels = plate_runner.raw_labels
        for i, label in enumerate(labels):
            if label in my_ignore_labels:
                continue
            
            strain = self.strains.setdefault(label, StrainData(label))
            strain.AddConditionRawCultureLevels(
                condition, plate_runner.raw_culture_levels[i,:])
            strain.AddConditionRawReporterLevels(
                condition, plate_runner.raw_reporter_levels[i,:])
            strain.AddConditionRawTimes(
                condition, plate_runner.raw_times[i,:])
            
        # Add filtered/computed data
        kept_labels = plate_runner.filtered_labels
        for i, label in enumerate(kept_labels):
            if label in my_ignore_labels:
                continue
            
            strain = self.strains.setdefault(label, StrainData(label))
            strain.AddConditionCultureLevels(
                condition, plate_runner.filtered_culture_levels[i,:])
            strain.AddConditionReporterLevels(
                condition, plate_runner.filtered_reporter_levels_no_bg[i,:])
            strain.AddConditionTimes(
                condition, plate_runner.shifted_times[i, :])
            
            activities = plate_runner.filtered_activities[i, :]
            smooth_activities = plate_runner.smooth_filtered_activities[i, :]
            max_activity = plate_runner.filtered_max_activities[i]
            strain.AddConditionActivities(condition, activities, max_activity)
            strain.AddConditionSmoothActivities(condition, smooth_activities, max_activity)
    
    def MakeStrainFigures(self, dirname):
        background_strain = self.strains[self.background_label]
        for s in self.strains.itervalues():
            s.MakeFigures(dirname, background_strain)

    def GetMeanMaxActivities(self, labels, condition):
        activities = []
        intervals = []
        for l in labels:
            s = self.strains[l]
            ac, err = s.GetMeanMaxActivity(condition)
            activities.append(ac)
            intervals.append(err)
        return activities, intervals
    
    def MakePerPlateFigures(self, dirname):
        fnames = []
        for condition, plates in self.plates.iteritems():
            for plate in plates:
                
                labels = plate.filtered_labels
                order = numpy.argsort(labels)
                smooth_activity = plate.smooth_filtered_activities
                max_activity = plate.filtered_max_activities
                left_mat = numpy.diag(1/max_activity)
                scaled_activity = numpy.dot(left_mat, smooth_activity)

                scaled_culture = plate.scaled_culture_levels
                max_culture_idx = numpy.argmin(numpy.abs(scaled_culture - 1.0), 1)
                for i, j in enumerate(max_culture_idx):
                    scaled_culture[i,j+1:] = 100
                
                Nr, _ = scaled_activity.shape
                od_fracs = numpy.arange(0.0, 1.0, 0.001)
                activity_per_od_frac = numpy.zeros((Nr, len(od_fracs)))
                for j, frac in enumerate(od_fracs):
                    abs_min = numpy.abs(scaled_culture - frac)
                    idxs = numpy.argmin(abs_min, 1)
                    for i in order:
                        idx = idxs[i]
                        activity_per_od_frac[i,j] = scaled_activity[i, idx]
                
                pylab.figure()
                pylab.title(condition)
                pylab.imshow(activity_per_od_frac, aspect='auto')
                
                condition_plate_name = util.slugify(condition)
                fname = '%s.png' % condition_plate_name
                pylab.savefig(path.join(dirname, fname),
                              format='png')
                
                fnames.append(fname)
            
        return fnames
        
    
    def MakeSummaryFigures(self, dirname, condition1, condition2):
        labels = self.GetStrainLabels()
        condition1_maxes, condition1_errs = self.GetMeanMaxActivities(
            labels, condition1)
        condition2_maxes, condition2_errs = self.GetMeanMaxActivities(
            labels, condition2)
        
        # Summary of all points
        pylab.figure()
        pylab.title('Maximal Activity of All Measured Strains')
        pylab.xlabel('Maximum in %s' % condition1)
        pylab.ylabel('Maximum in %s' % condition2)
        pylab.loglog(condition1_maxes, condition2_maxes, 'b.', basey=2)
        pylab.errorbar(condition1_maxes, condition2_maxes, yerr=condition2_errs,
                       xerr=condition1_errs, fmt=None)
        
        log_x = numpy.log2(condition1_maxes)
        log_y = numpy.log2(condition2_maxes)
        min_x, max_x = numpy.min(log_x), numpy.max(log_x)
        log_x_range = numpy.arange(min_x, max_x, 0.5)
        y_int = st.nanmean(log_y) - st.nanmean(log_x)
        predicted_y = log_x + y_int
        diffs = log_y - predicted_y
        res = abs(diffs)
        r2 = util.calc_r2(log_y, predicted_y)
        pylab.loglog(numpy.exp2(log_x), numpy.exp2(log_x), 'k--',
                     basey=2, basex=2,
                     label='y = x')
        pylab.loglog(numpy.exp2(log_x_range), numpy.exp2(log_x_range + y_int), 'b--',
                     basey=2, basex=2,
                     label='y = x + %.2g (r2 = %.2g)' % (y_int, r2))
        
        # ~28% induction is the minimum.
        far_points = pylab.find(res > 0.4)
        #for i in far_points:
        for i in xrange(len(res)):
            color = GetGeneTextColor(labels[i])
            pylab.text(condition1_maxes[i], condition2_maxes[i],
                       labels[i], color=color)
        
        pylab.legend(loc='lower right', prop=LEGEND_FONT)
        summary_fname = 'summary.png'
        pylab.savefig(path.join(dirname, summary_fname),
                      format='png')
        
        # Summary of most differentially expressed points
        pylab.figure()
        pylab.title('Maximal Activity of Differentially Expressed Strains')
        pylab.xlabel('Maximum in %s' % condition1)
        pylab.ylabel('Maximum in %s' % condition2)
        
        far_logx = log_x[far_points]
        far_logy = log_y[far_points]
        big_diffs = far_logy - far_logx
        above_line = pylab.find(big_diffs > 0)
        below_line = pylab.find(big_diffs < 0)
        
        pylab.loglog(numpy.exp2(log_x_range), numpy.exp2(log_x_range), 'k--',
                     label='y = x', basey=2, basex=2)
        
        # Plot points lying above the line
        above_line_logx = far_logx[above_line]
        above_line_logy = far_logy[above_line]
        pylab.loglog(numpy.exp2(above_line_logx), numpy.exp2(above_line_logy),
                     'g.', basey=2, basex=2)
        
        y_int = pylab.mean(above_line_logy) - pylab.mean(above_line_logx)
        predicted_y = above_line_logx + y_int
        r2 = util.calc_r2(above_line_logy, predicted_y)
        pylab.loglog(numpy.exp2(log_x_range), numpy.exp2(log_x_range + y_int),
                     'g--', basey=2, basex=2,
                     label='f=x + %.2g (r2=%.2f)' % (y_int, r2))
        
        # Plot points lying below the line.
        below_line_logx = far_logx[below_line]
        below_line_logy = far_logy[below_line]
        pylab.loglog(numpy.exp2(below_line_logx), numpy.exp2(below_line_logy),
                     'r.', basey=2, basex=2)
        
        y_int = pylab.mean(below_line_logy) - pylab.mean(below_line_logx)
        predicted_y = below_line_logx + y_int
        r2 = util.calc_r2(below_line_logy, predicted_y)
        pylab.loglog(numpy.exp2(log_x_range), numpy.exp2(log_x_range + y_int),
                     'r--', basey=2, basex=2,
                     label='f=x + %.2g (r2=%.2f)' % (y_int, r2))
        
        for i in far_points:
            color = GetGeneTextColor(labels[i])
            pylab.text(condition1_maxes[i], condition2_maxes[i],
                       labels[i], color=color)
            
        pylab.legend(loc='lower right', prop=LEGEND_FONT)
        diff_fname = 'differential.png'
        pylab.savefig(path.join(dirname, diff_fname),
                      format='png')
        
        return [summary_fname, diff_fname]
        
        
def Main():
    options, _ = MakeOpts().parse_args(sys.argv)
    assert options.experiment_id
    assert options.first_plate_ids and options.second_plate_ids
    assert options.culture_label and options.reporter_label
    assert options.output_dir
    
    if not path.exists(options.output_dir):
        util._mkdir(options.output_dir)
    
    imgs_path = path.join(options.output_dir, 'imgs/')
    if not path.exists(imgs_path):
        util._mkdir(imgs_path)

    first_plate_ids = map(str.strip, options.first_plate_ids.split(','))
    second_plate_ids = map(str.strip, options.second_plate_ids.split(','))
    
    labels_to_ignore = set()
    for l in options.labels_to_ignore.split(','):
        labels_to_ignore.add(l.strip())
    
    print 'Reading plates from experiment %s' % (options.experiment_id)
    db = MySQLDatabase(host='hldbv02', user='ronm', 
                       passwd='a1a1a1', db='tecan')

    filterer = promoter_activity.CultureReporterFilterer(options.min_culture_level,
                                                         options.min_reporter_level)
    reporter_bg_sub = promoter_activity.ReporterBackgroundSubtracter(
        options.background_label)
    culture_shifter = promoter_activity.CultureShifter()
    activity_calc = promoter_activity.ReporterActivityCalculator(
        options.lower_culture_bound, options.upper_culture_bound,
        min_reporter_level=options.lower_reporter_bound,
        window_size=options.window_size)

    first_plate_runners = []
    second_plate_runners = []
    print 'Calculating promoter activities for first condition'
    runner1 = PlateActivityRunner(
        options.culture_label, options.reporter_label,
        filterer, culture_shifter, reporter_bg_sub, activity_calc)
    
    for plate_id in first_plate_ids:
        plate = Plate96.FromDatabase(db, options.experiment_id, plate_id)
        runner1.AddPlate(plate)
    
    runner1.Run()
    first_plate_runners.append(runner1)

    print 'Calculating promoter activities for second condition'
    runner2 = PlateActivityRunner(
        options.culture_label, options.reporter_label,
        filterer, culture_shifter, reporter_bg_sub, activity_calc)
    
    for plate_id in second_plate_ids:
        plate = Plate96.FromDatabase(db, options.experiment_id, plate_id)
        runner2.AddPlate(plate)
    
    runner2.Run()
    second_plate_runners.append(runner2)
    
    # Unify strain data.
    print 'Saving figures'
    strains_data = StrainConditionsData(options.background_label)
    for plate_data in first_plate_runners:
        strains_data.AddPlateData('Glucose', plate_data,
                                  ignore_labels=labels_to_ignore)
    for plate_data in second_plate_runners:
        strains_data.AddPlateData('Gluconate', plate_data,
                                  ignore_labels=labels_to_ignore)
    strains_data.MakeStrainFigures(imgs_path)
    summary_fignames = strains_data.MakeSummaryFigures(
        imgs_path, 'Glucose', 'Gluconate')
    plate_fignames = strains_data.MakePerPlateFigures(imgs_path)
    
    
    labels = strains_data.GetStrainLabels()
    condition1_activities, condition1_errs = strains_data.GetMeanMaxActivities(
        labels, 'Glucose')
    condition2_activities, condition2_errs = strains_data.GetMeanMaxActivities(
        labels, 'Gluconate')
    log_1 = numpy.log2(condition1_activities)
    log_2 = numpy.log2(condition2_activities)
    diffs = log_2 - log_1
    sorted_diffs = list(numpy.argsort(diffs))
    sorted_diffs.reverse()
    diffs_data = []
    for i in sorted_diffs:
        logfold = diffs[i]
        fold = numpy.exp2(logfold)
        if numpy.isnan(logfold):
            logfold = None
            fold = None
        
        diffs_data.append({'label': labels[i],
                           'fold_change': fold,
                           'log_fold': logfold})        
    
    # Render the template.
    print 'Writing HTML output'
    template_data = {'experiment_id': options.experiment_id,
                     'first_plate_ids': first_plate_ids,
                     'second_plate_ids': second_plate_ids,
                     'culture_label': options.culture_label,
                     'reporter_label': options.reporter_label,
                     'first_plates': first_plate_runners,
                     'second_plates': second_plate_runners,
                     'strains_data': strains_data,
                     'diffs_data': diffs_data,
                     'summary_figure_fnames': summary_fignames,
                     'per_plate_figure_fnames': plate_fignames}
    template_fname = path.join(options.output_dir, 'results.html')
    templates.render_to_file(
        'compare_promoter_activities.html', template_data, template_fname)
    
    return
        
    
    
if __name__ == '__main__':
    Main()