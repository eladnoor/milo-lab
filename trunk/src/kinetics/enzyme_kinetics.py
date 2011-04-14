#!/usr/bin/python

import logging 
import pylab

from scipy.stats import stats
from toolbox import database


class EnzymeKineticsData(object):
    
    def __init__(self, ec,
                 reaction_id,
                 substrate_id,
                 molecular_weight,
                 specific_activity):
        self.ec = ec
        self.rid = reaction_id
        self.cid = substrate_id
        self.mw = molecular_weight
        self.sa = specific_activity

    @staticmethod
    def FromDBRow(row):
        """Initialize from a database row."""
        ec = row['ec']
        rid = row['rid']
        cid = row['cid']
        mw = row['mw']
        sa = row['sa']
        return EnzymeKineticsData(ec, rid, cid, mw, sa)
    

class EnzymeKineticsContainer(object):
    """Contains multiple EnzymeKineticsData per EC."""
    
    def __init__(self, ec,
                 data=None):
        self.data = data or []
        self.ec = ec
    
    def AddEnzymeKineticsData(self, data_obj):
        """Add a data to the collection.
        
        Returns a self-reference for fun.
        """
        self.data.append(data_obj)
        return self
    
    def GetKineticsForSubstrate(self, cid):
        """Returns an EnzymeKineticsContainer for the given substrate.
        
        Filters out all the measurements that are not for the given 
        substrate.
        
        Args:
            cid: the KEGG ID of the substrate.
        """
        matches_cid = lambda x: x.cid == cid
        matching_data = filter(matches_cid, self.data)
        return EnzymeKineticsContainer(self.ec, matching_data)
    
    def AllSpecificActivities(self):
        """Returns a list of all specific activities contained."""
        get_sa = lambda x: x.sa
        return map(get_sa, self.data)        
    
    def MeanSpecificActivity(self):
        """Get the mean specific activity of all kinetic measurements."""
        return pylab.mean(self.AllSpecificActivities())
    
    def HarmonicMeanSpecificActivity(self):
        """Get the mean specific activity of all kinetic measurements."""
        return stats.hmean(self.AllSpecificActivities())
    
    @staticmethod
    def FromDBRow(row):
        """Initialize from a database row.
        
        Creates the container and adds the row.
        """
        ec = row['ec']
        datum = EnzymeKineticsData.FromDBRow(row)
        return EnzymeKineticsContainer(ec).AddEnzymeKineticsData(datum)
    
    
class EnzymeKinetics(object):
    
    def __init__(self, db):
        self.db = db
        self._kinetics_by_ec = {}
        self._LoadFromDB()
    
    def _LoadFromDB(self):
        """Loads data from the database."""
        for i, row in enumerate(self.db.DictReader('merged_sa')):
            ec = row['ec']
            if not ec:
                logging.warning('Missing EC # for row %s', row)
                continue
            
            if i % 25 == 0:
                logging.info('Adding entry for EC %s', ec)

            container = self._kinetics_by_ec.get(ec)
            if container:
                container.AddEnzymeKineticsData(EnzymeKineticsData.FromDBRow(row))
            else:
                container = EnzymeKineticsContainer.FromDBRow(row)
                self._kinetics_by_ec[ec] = container
    
    def GetKinetics(self, ec):
        """Returns kinetic information for this EC number."""
        return self._kinetics_by_ec.get(ec)


if __name__ == '__main__':
    db = database.SqliteDatabase('/home/flamholz/Dropbox/BRENDA/enzyme_db.sqlite')
    kinetics = EnzymeKinetics(db)
    print 'Kinetics!'
    hexokinase = kinetics.GetKinetics('2.7.1.1')
    hexo_glucose = hexokinase.GetKineticsForSubstrate(31)
    print 'Harmonic Mean %f' % hexo_glucose.HarmonicMeanSpecificActivity()
