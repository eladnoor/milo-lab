#!/usr/bin/python

"""Functions relating to our different thermodynamic_esimators."""

import logging

from toolbox.database import SqliteDatabase
from pygibbs.groups import GroupContribution
from pygibbs.hatzimanikatis import Hatzi
from pygibbs.thermodynamics import PsuedoisomerTableThermodynamics
from pygibbs.thermodynamics import BinaryThermodynamics
from pygibbs.thermodynamics import ReactionThermodynamics

ESTIMATOR_NAMES = ('hatzi_gc', 'BGC', 'PGC', 'merged')

def EstimatorNames():
    return ESTIMATOR_NAMES


def LoadAllEstimators():
    db_public = SqliteDatabase('../data/public_data.sqlite')
    db_gibbs = SqliteDatabase('../res/gibbs.sqlite')
    tables = {'alberty': (db_public, 'alberty_pseudoisomers', 'Alberty'),
              'PRC': (db_gibbs, 'prc_pseudoisomers', 'our method (PRC)')}

    estimators = {}

    for key, (db, table_name, thermo_name) in tables.iteritems():
        if db.DoesTableExist(table_name):
            estimators[key] = PsuedoisomerTableThermodynamics.FromDatabase(
                                            db, table_name, name=thermo_name)
        else:
            logging.warning('The table %s does not exist in %s' % (table_name, str(db)))
    
    estimators['hatzi_gc'] = Hatzi(use_pKa=False)
    #estimators['hatzi_gc_pka'] = Hatzi(use_pKa=True)
    
    estimators['BGC'] = GroupContribution(db=db_gibbs, transformed=True)
    estimators['BGC'].init()
    estimators['BGC'].name = 'our method (BGC)'

    estimators['PGC'] = GroupContribution(db=db_gibbs, transformed=False)
    estimators['PGC'].init()
    estimators['PGC'].name = 'our method (PGC)'
    
    estimators['merged'] = BinaryThermodynamics(estimators['alberty'],
                                                estimators['PGC'])
    
    estimators['C1'] = ReactionThermodynamics.FromCsv(
        '../data/thermodynamics/c1_reaction_thermodynamics.csv',
        estimators['alberty'])
    
    estimators['merged_C1'] = BinaryThermodynamics(estimators['C1'],
                                                   estimators['PGC'])

    for thermo in estimators.values():
        thermo.load_bounds('../data/thermodynamics/concentration_bounds.csv')

    return estimators