import logging
import pylab
import re

from django.db import models
from group_contribution import constants


class CommonName(models.Model):
    name = models.CharField(max_length=500)
    enabled = models.BooleanField()
    
    def __unicode__(self):
        return self.name
    

class Specie(models.Model):
    """A single specie of a compound."""
    # The ID of the compound in KEGG.
    kegg_id = models.CharField(max_length=10)
    
    # The number of hydrogens in the species.
    number_of_hydrogens = models.IntegerField()
    
    # The net charge (eV).
    net_charge = models.IntegerField()


class ValueSource(models.Model):
    """The source of a particular numeric value."""
    # The name of the source.
    name = models.CharField(max_length=100)
    
    # A link explaining the source.
    link = models.URLField(null=True)
    
    
class SpeciesFormationEnergy(models.Model):
    """The standard formation energy of a single species."""
    # Which specie this value is for.
    specie = models.ForeignKey(Specie)
    
    # The standard formation energy in kJ/mol.
    value = models.FloatField()
    
    # The source of this value.
    source = models.ForeignKey(ValueSource)
    
    def transform(self,
                  pH=constants.DEFAULT_PH,
                  ionic_strength=constants.DEFAULT_IONIC_STRENGTH,
                  temp=constants.DEFAULT_TEMP):
        """Transform this individual estimate to difference conditions."""
        # Short names are nice!
        _i_s = ionic_strength
        _r = constants.R
        _n_h = self.specie.number_of_hydrogens
        _n_c = self.specie.net_charge
        
        chem_potential = _n_h * _r * temp * pylab.log(10) * pH
        ionic_potential = (2.91482 * (_n_c ** 2 - _n_h) * pylab.sqrt(_i_s) /
                           (1 + 1.6 * pylab.sqrt(_i_s)))
        return self.value + chem_potential - ionic_potential


class Compound(models.Model):
    # The ID of the compound in KEGG.
    kegg_id = models.CharField(max_length=10)
    
    # A list of common names of the compound, used for searching.
    common_names = models.ManyToManyField(CommonName)
    
    # The chemical formula.
    formula = models.CharField(max_length=500)
    
    # The molecular mass.
    mass = models.FloatField()  # In Daltons.
    
    # Estimates of Delta G for this compound.
    species_formation_energies = models.ManyToManyField(SpeciesFormationEnergy)

    def GetFormationEnergy(self, pH=constants.DEFAULT_PH,
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
        all_sfes = self.species_formation_energies.all()
        if not all_sfes:
            # No data...
            return None
        
        # Shorter names are handy!
        _i_s = ionic_strength
        _r = constants.R
        
        # Compute per-species transforms, scaled down by R*T.
        transform = lambda x: x.transform(pH, _i_s, temp)
        scaled_transforms = [(-transform(sfe) / (_r * temp)) for sfe in all_sfes]
        
        # Numerical issues: taking a sum of exp(v) for |v| quite large.
        # Use the fact that we take a log later to offset all values by a 
        # constant (the minimum value).
        offset = min(scaled_transforms)
        scaled_offset_transforms = [(st - offset) for st in scaled_transforms]
        sum_exp = sum(pylab.exp(scaled_offset_transforms))
        return - _r * temp * (offset + pylab.log(sum_exp))

    def GetAtomBag(self):
        """Returns a dictionary of atoms and their counts for this compound."""
        if not self.formula:
            logging.error('Formula is not defined for KEGG ID %s', self.kegg_id)
            return None
        
        if '(' in self.formula or ')' in self.formula:
            return None
            
        atom_bag = {}
        for atom, count in re.findall("([A-Z][a-z]*)([0-9]*)", self.formula):
            count = count or 1
            count = int(count)
            atom_bag[atom] = count
        
        # Wildcards are not allowed.
        if 'R' in atom_bag:
            return None
        
        return atom_bag

    def GetKeggLink(self):
        """Returns a link to the KEGG page for this compound."""
        if not self.kegg_id:
            return None
        
        return 'http://kegg.jp/dbget-bin/www_bget?cpd:%s' % self.kegg_id

    def GetHtmlFormattedFormula(self):
        """Returns the chemical formula with HTML formatted subscripts."""
        if not self.formula:
            return None
        
        return re.sub(r'(\d+)', r'<sub>\1</sub>', self.formula)
    
    html_formula = property(lambda self: self.GetHtmlFormattedFormula())
    kegg_link = property(lambda self: self.GetKeggLink())
    all_common_names = property(lambda self: self.common_names.all())
    all_formation_energies = property(
        lambda self: self.species_formation_energies.all())
    
    
    def __unicode__(self):
        names = self.common_names.all()
        if names:
            return unicode(names[0])
        return unicode(self.formula)
    
    @staticmethod
    def _GetAtomDiff(collection):
        """Get the net atom counts from the collection.
        
        Args:
            collection: an iterable of 2 tuples (coeff, kegg_id).
                Coefficients can be negative.
        """
        ids = [x[1] for x in collection]
        compounds = Compound.GetCompoundsByKeggId(ids)

        atom_diff = {}
        for coeff, id in collection:
            if id not in compounds:
                logging.error('Unknown compound id %s', id)
                return None
            
            c = compounds[id]
            atom_bag = c.GetAtomBag()
            if not atom_bag:
                logging.error('Failed to fetch atom bag for %s', c.formula)
                return None
            
            for atomic_number, atom_count in atom_bag.iteritems():
                new_diff = atom_diff.get(atomic_number, 0) + coeff * atom_count
                atom_diff[atomic_number] = new_diff
        
        return atom_diff
        
    @staticmethod
    def _IsBalanced(atom_diff):
        """Checks if the per-atom diffs represent a balanced collection.
        
        Args:
            atom_diff: a dictionary mapping atomic numbers to counts.
            
        Returns:
            True if balanced.
        """
        if not atom_diff:
            return False
        
        return max([abs(x) for x in atom_diff.values()]) < 0.01
    
    @staticmethod
    def IsBalanced(collection):
        """Checks if the collection is atom-wise balanced.
        
        Args:
            collection: an iterable of 2 tuples (coeff, kegg_id).
                Coefficients can be negative.
        
        Returns:
            True if the collection is atom-wise balanced.
        """
        return Compound._IsBalanced(Compound._GetAtomDiff(collection))        
    
    @staticmethod
    def TryBalanceWithWater(collection):
        """Tries to balance the reaction with water only.
        
        Args:
            Args: an iterable of 2 tuples (coeff, kegg_id).
                Coefficients can be negative.
        
        Returns:
            The number of waters to add/subtract to the reaction or None if impossible.
        """
        atom_diff = Compound._GetAtomDiff(collection)
        if not atom_diff:
            return None
        
        # Ignore hydrogen.
        atom_diff.pop('H')
        
        # Omit oxygen for checking balancedness.
        oxy_count = atom_diff.pop('O', 0)
        
        # If it's still not balanced, then there's trouble.
        if not Compound._IsBalanced(atom_diff):
            return None
        
        # If it's balanced without considering oxygen/hydrogen,
        # then we need to add/subtract whatever the oxygen excess is.
        return (- oxy_count)
        
    @staticmethod
    def ReactionIsBalanced(reactants, products):
        """Returns True if the reaction is balanced.
        
        Args:
            reactants: an iterable of 2 tuples (coeff, kegg_id) for reactants.
            products: an iterable of 2 tuples (coeff, kegg_id) for products.
        """
        collection = reactants + [(-coeff, id) for coeff, id in products]
        return Compound.IsBalanced(collection)
    
    @staticmethod
    def _AddWaterTo(reaction_side, water_amt):
        """Adds an amount of water to a single side of a reaction.
        
        Args:
            reaction_side: An iterable of 2 tuples (coeff, kegg_id).
            water_amt: the amount of water to add.
        
        Returns:
            A list of 2 tuples with water_amt added.
        """
        water_id = 'C00001'
        side = list(reaction_side)
        found_water = False
        for i, tpl in enumerate(side):
            coeff, id = tpl
            if id == water_id:
                found_water = True
                coeff += water_amt
                side[i] = (coeff, id)
                
        if not found_water:
            side.append((water_id, water_amt))
            
        return side                
    
    @staticmethod
    def CanBalanceReactionWithWater(reactants, products):
        """Attempts to balance the reaction with water.
        
        Args:
            reactants: an iterable of 2 tuples (coeff, kegg_id) for reactants.
            products: an iterable of 2 tuples (coeff, kegg_id) for products.
            
        Returns:
            True if the reaction can be balanced with water.
        """
        collection = reactants + [(-coeff, id) for coeff, id in products]
        water_excess = Compound.TryBalanceWithWater(collection)
        if water_excess == None:
            return False
        return True
        
    @staticmethod
    def GetTotalFormationEnergy(collection,
                                pH=constants.DEFAULT_PH,
                                ionic_strength=constants.DEFAULT_IONIC_STRENGTH,
                                temp=constants.DEFAULT_TEMP):
        """Compute an estimate for a collection of compounds + coefficients.
        
        You can compute the DeltaG of a reaction using negative coefficients for
        products.
        
        Args:
            collection: an iterable of 2 tuples (coeff, kegg_id).
                Coefficients can be negative. 
        """
        ids = [x[1] for x in collection]
        compounds = Compound.GetCompoundsByKeggId(ids)
        
        sum = 0
        for coeff, id in collection:
            if id not in compounds:
                logging.info('Could not find compound %s', id)
                return None
            
            c = compounds[id]
            est = c.GetFormationEnergy(pH, ionic_strength, temp)
            if not est:
                logging.info('No estimate for compound %s', id)
                return None
            
            sum += coeff * est
        
        return sum
    
    @staticmethod
    def GetReactionEnergy(reactants, products,
                          pH=constants.DEFAULT_PH,
                          ionic_strength=constants.DEFAULT_IONIC_STRENGTH,
                          temp=constants.DEFAULT_TEMP):
        """Compute the DeltaG for a reaction.
        
        Args:
            reactants: an iterable of 2 tuples (coeff, kegg_id) for reactants.
            products: an iterable of 2 tuples (coeff, kegg_id) for products.
        
        Returns:
            The DeltaG for the reaction, or None if data was missing.
        """
        reactants_sum = Compound.GetTotalFormationEnergy(
            reactants, pH, ionic_strength, temp)
        products_sum = Compound.GetTotalFormationEnergy(
            products, pH, ionic_strength, temp)
        if not products_sum:
            logging.warning('Failed to get products formation energy.')
            return None
        if not reactants_sum:
            logging.warning('Failed to get reactants formation energy.')
            return None
        
        return products_sum - reactants_sum
    
    @staticmethod
    def AsLibrary():
        """Get's a map of common names to Compound objects from the DB."""
        all_compounds = Compound.objects.all()
        names_to_compounds = {}
        for compound in all_compounds:
            for name in compound.common_names.all():
                names_to_compounds[name] = compound
        
        return names_to_compounds
    
    @staticmethod
    def GetCompoundsByKeggId(kegg_ids):
        """Fetch compounds from a list of KEGG IDs.
        
        Args:
            kegg_ids: a list of KEGG IDs.
        
        Returns:
            A dictionary mapping KEGG ID to Compounds.
        """
        compounds = Compound.objects.filter(kegg_id__in=kegg_ids)
        return dict((c.kegg_id, c) for c in compounds if c != None)
        
        
    
