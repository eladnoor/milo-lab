#!/usr/bin/python

import constants
import logging
import pylab

##########################################################################
## TODO(flamholz): All classes in this file should be backed by a database.
##########################################################################

class CompoundEstimate(object):
    """An estimate for a particular compound."""
    
    def __init__(self, kegg_id, protonation_level, net_charge, dG0_estimate):
        """Construct a CompoundEstimate from data.
        
        Args:
            kegg_id: the KEGG id of the compound.
            protonation_level: the protonation level for this estimate.
            net_charge: the net charge of this protonation level.
            dG0_estimate: the delta G_0 estimate for this protonation level/net charge.
        """
        self.kegg_id = kegg_id
        self.protonation_level = protonation_level
        self.net_charge = net_charge
        self.dG0_estimate = dG0_estimate
    
    def transform(self,
                  pH=constants.DEFAULT_PH,
                  ionic_strength=constants.DEFAULT_IONIC_STRENGTH,
                  temp=constants.DEFAULT_TEMP):
        chem_potential = self.protonation_level * constants.R * temp * pylab.log(10) * pH
        ionic_potential = (2.91482 * (self.net_charge ** 2 - self.protonation_level) * pylab.sqrt(ionic_strength) /
                           (1 + 1.6 * pylab.sqrt(ionic_strength)))
        return self.dG0_estimate + chem_potential - ionic_potential
    
    @staticmethod
    def FromCsvLine(line):
        """Reads a compound estimate from my CSV file.
    
        Arguments:
            line: the line of CSV to read from.
        """
        fields = line.strip().split(', ')
        if len(fields) != 4:
            logging.error('Failed to parse the CSV line.')
            return None
        
        try:
            kegg_id = 'C%05d' % int(fields[0])
            protonation_level = int(fields[1])
            net_charge = int(fields[2])
            delta_g0_est = float(fields[3]) 
            return CompoundEstimate(kegg_id, protonation_level,
                                    net_charge, delta_g0_est)
        except ValueError, e:
            logging.error(e)
        
        return None
    

class CompoundEstimates(object):
    """A collection of CompoundEstimates."""
    
    def __init__(self):
        self._estimates = {}
    
    def AddEstimate(self, compound_estimate):
        """Adds a single estimate to the collection.
        
        Args:
            compound_estimate: the CompoundEstimate object to add.
        """
        if not compound_estimate:
            return
        self._estimates.setdefault(compound_estimate.kegg_id,
                                   []).append(compound_estimate)
    
    def GetCompoundEstimate(self,
                            kegg_id,
                            pH=constants.DEFAULT_PH,
                            ionic_strength=constants.DEFAULT_IONIC_STRENGTH,
                            temp=constants.DEFAULT_TEMP):
        """Get a deltaG estimate for the given compound.
        
        Args:
            kegg_id: the KEGG id of the compound.
            pH: the PH to estimate at.
            ionic_strength: the ionic strength to estimate at.
            temp: the temperature to estimate at.
        
        Returns:
            The estimated delta G in the given conditions or None.
        """
        if kegg_id not in self._estimates:
            return None
        
        # Scale transforms down by R*T.
        scaled_transforms = [(ce.transform() / (constants.R * temp))
                             for ce in self._estimates[kegg_id]]
        
        # Numerical issues: taking a sum of exp(v) for |v| quite large.
        # Use the fact that we take a log later to offset all values by a 
        # constant (the minimum value).
        offset = min(scaled_transforms)
        scaled_offset_transforms = [(st - offset) for st in scaled_transforms]
        sum_exp = sum(pylab.exp(scaled_offset_transforms))
        return constants.R * temp * (offset + pylab.log(sum_exp))
    
    def _GetCollectionEstimate(self, collection,
                               pH=constants.DEFAULT_PH,
                               ionic_strength=constants.DEFAULT_IONIC_STRENGTH,
                               temp=constants.DEFAULT_TEMP):
        """Compute an estimate for a collection of compounds + coefficients.
        
        Args:
            collection: an iterable of 2 tuples (coeff, kegg_id).
        """
        sum = 0.0
        for coeff, compound_id in collection:
            est = self.GetCompoundEstimate(compound_id, pH,
                                           ionic_strength, temp)
            if est == None:
                return None
            
            sum += coeff * est
            
        return sum
    
    def GetReactionEstimate(self, reactants, products,
                            pH=constants.DEFAULT_PH,
                            ionic_strength=constants.DEFAULT_IONIC_STRENGTH,
                            temp=constants.DEFAULT_TEMP):
        """Compute an estimate for a reaction.
        
        Args:
            reactants: an iterable of 2 tuples (coeff, kegg_id) for reactants.
            products: an iterable of 2 tuples (coeff, kegg_id) for products.
        """
        reactants_sum = self._GetCollectionEstimate(reactants, pH, ionic_strength, temp)
        products_sum = self._GetCollectionEstimate(products, pH, ionic_strength, temp)
        if not products_sum or not reactants_sum:
            return None
        
        return products_sum - reactants_sum
    
    @staticmethod
    def FromCsvFile(filename):
        """Build a collection for CompoundEstimates from a CSV file.
        
        Args:
            filename: the name of the csv file.
        """
        ces = CompoundEstimates()
        f = open(filename, 'r')
        for line in f:
            # Hack: ignore header line.
            if line.startswith('cid'):
                continue
            
            c = CompoundEstimate.FromCsvLine(line)
            ces.AddEstimate(c)
        
        f.close()
        return ces

if __name__ == '__main__':
    ces = CompoundEstimates.FromCsvFile('kegg_thermo_data.csv')
    print "ATP:", ces.GetCompoundEstimate('C00002')
    print "Glucosamine:", ces.GetCompoundEstimate('C04886')