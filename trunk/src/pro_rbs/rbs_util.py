#!/usr/bin/python

import logging
import pylab

from pro_rbs.rbs_calc.RBS_Calculator import RBS_Calculator

R = 8.31e-3 # kJ/(K*mol)
JOULES_PER_CAL = 4.184
TEMP = 298.15 # K


def TryCalcDG(seq):
    """Wrapper around weird library to calculate binding dG."""
    seq_str = str(seq)
    start_range = [0, len(seq_str)]
    name = seq_str
    calc_obj = RBS_Calculator(seq_str, start_range, name)
    calc_obj.calc_dG()
    if not calc_obj.dG_total_list:
        logging.warning('Failed to calculate dG')
        raise Exception('Failed to calculate dG')

    return calc_obj.dG_total_list[0]    
    

def TokJ(dG):
    return JOULES_PER_CAL * dG


def PredictedRate(dG):
    dG_kj = JOULES_PER_CAL * dG
    rt = R*TEMP
    return pylab.exp(-dG/rt)


def LogRate(dG):
    return pylab.log10(PredictedRate(dG))