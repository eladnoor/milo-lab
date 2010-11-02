import logging
import pylab
import re

from django.db import models
from gibbs import constants


class CommonName(models.Model):
    """A common name of a compound."""
    name = models.CharField(max_length=500)
    enabled = models.BooleanField()
    
    def __unicode__(self):
        return self.name
    
    
class ValueSource(models.Model):
    """The source of a particular numeric value."""
    # The name of the source.
    name = models.CharField(max_length=100)
    
    # A link explaining the source.
    link = models.URLField(null=True)
    
    def __unicode__(self):
        return self.name


class SpeciesFormationEnergy(models.Model):
    """The standard formation energy of a single species."""
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

    def __unicode__(self):
        return str(self.value)
    

class Specie(models.Model):
    """A single specie of a compound."""
    # The ID of the compound in KEGG.
    kegg_id = models.CharField(max_length=10)
    
    # The number of hydrogens in the species.
    number_of_hydrogens = models.IntegerField()
    
    # The net charge (eV).
    net_charge = models.IntegerField()
    
    # The formation energy of this specie, measured or approximated.
    formation_energy = models.OneToOneField(SpeciesFormationEnergy)
    
    def __unicode__(self):
        return self.kegg_id


class Compound(models.Model):
    """A single compound."""
    # The ID of the compound in KEGG.
    kegg_id = models.CharField(max_length=10, null=True)
    
    # A list of common names of the compound, used for searching.
    common_names = models.ManyToManyField(CommonName)
    
    # The chemical formula.
    formula = models.CharField(max_length=500)
    
    # The molecular mass.
    mass = models.FloatField()  # In Daltons.
    
    # Estimates of Delta G for this compound.
    species = models.ManyToManyField(Specie)
    
    # An explanation for when no DeltaG estimate is available.
    no_dg_explanation = models.CharField(max_length=2048,
                                         blank=True,
                                         null=True)

    def __init__(self, *args, **kwargs):
        super(Compound, self).__init__(*args, **kwargs)
        self._all_species = None
        
    def DeltaG(self, pH=constants.DEFAULT_PH,
               ionic_strength=constants.DEFAULT_IONIC_STRENGTH):
        """Get a deltaG estimate for the given compound.
        
        Args:
            pH: the PH to estimate at.
            ionic_strength: the ionic strength to estimate at.
            temp: the temperature to estimate at.
        
        Returns:
            The estimated delta G in the given conditions or None.
        """
        if not self.all_species:
            # No data...
            return None
        
        # Shorter names are handy!
        _i_s = ionic_strength
        _r = constants.R
        _t = constants.DEFAULT_TEMP

        # Compute per-species transforms, scaled down by R*T.
        transform = lambda x: x.transform(pH, _i_s, _t)
        scaled_transforms = [(-transform(s.formation_energy) / (_r * _t))
                             for s in self.all_species]
        
        # Numerical issues: taking a sum of exp(v) for |v| quite large.
        # Use the fact that we take a log later to offset all values by a 
        # constant (the minimum value).
        offset = min(scaled_transforms)
        scaled_offset_transforms = [(st - offset) for st in scaled_transforms]
        sum_exp = sum(pylab.exp(scaled_offset_transforms))
        return - _r * _t * (offset + pylab.log(sum_exp))

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

    def GetLink(self):
        if not self.kegg_id:
            return None
        return '/compound?compoundId=%s' % self.kegg_id

    def GetKeggLink(self):
        """Returns a link to the KEGG page for this compound."""
        if not self.kegg_id:
            return None
        
        return 'http://kegg.jp/dbget-bin/www_bget?cpd:%s' % self.kegg_id
    
    def GetSmallImageUrl(self):
        """Returns the URL of a small image of the compound structure."""
        if not self.kegg_id:
            return None
        
        return 'http://kegg.jp/Fig/compound_small/%s.gif' % self.kegg_id

    def GetHtmlFormattedFormula(self):
        """Returns the chemical formula with HTML formatted subscripts."""
        if not self.formula:
            return None
        
        return re.sub(r'(\d+)', r'<sub>\1</sub>', self.formula)
    
    def GetSpecies(self):
        """Gets the list of SpeciesFormationEnergies, potentially caching."""
        if not self._all_species:
            self._all_species = self.species.all()
        return self._all_species
    
    html_formula = property(GetHtmlFormattedFormula)
    link = property(GetLink)
    kegg_link = property(GetKeggLink)
    small_image_url = property(GetSmallImageUrl)
    all_common_names = property(lambda self: self.common_names.all())
    all_species = property(GetSpecies)
    standard_formation_energy = property(DeltaG)
    
    def StashTransformedSpeciesEnergies(self, ph, ionic_strength):
        """Stash the transformed species formation energy in each one."""
        for species in self.all_species:
            species.transformed_energy = species.formation_energy.transform(
                ph, ionic_strength)
    
    def __unicode__(self):
        """Return a single string identifier of this Compound."""
        names = self.all_common_names
        if names:
            return unicode(names[0])
        return unicode(self.formula)
    
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


# TODO(flamholz): The two classes below don't belong in this file...


class CompoundWithCoeff(object):
    """A compound with a stoichiometric coefficient."""
    
    def __init__(self, coeff, compound, name=None):
        self.compound = compound
        self.coeff = coeff
        self.name = name
        
    def Minus(self):
        """Returns a new CompoundWithCoeff with coeff = -self.coeff."""
        return CompoundWithCoeff(-self.coeff, self.compound, self.name)
    

class Reaction(models.Model):
    """A reaction."""
    
    def __init__(self, reactants=None, products=None):
        """Construction.
        
        Args:
            reactants: a list of CompoundWithCoeff objects.
            products: a list of CompoundWithCoeff objects.
        """
        self.reactants = reactants or []
        self.products = products or []
    
    @staticmethod
    def FromIds(reactants, products):
        """Build a reaction object from lists of IDs.
        
        Args:
            reactants: an iterable of (coeff, kegg_id, name) of reactants.
            products: an iterable of (coeff, kegg_id, name) of products.
            
        Returns:
            A properly set-up Reaction object or None if there's an error.
        """
        r_ids = [id for unused_coeff, id, unused_name in reactants]
        p_ids = [id for unused_coeff, id, unused_name in products]
        compounds = Compound.GetCompoundsByKeggId(r_ids + p_ids)
        
        # Build the reaction object.
        rxn = Reaction()        
        for coeff, id, name in reactants:
            if id not in compounds:
                logging.error('Unknown reactant %s', id)
                return None
            
            rxn.reactants.append(CompoundWithCoeff(coeff, compounds[id], name))
        
        for coeff, id, name in products:
            if id not in compounds:
                logging.error('Unknown product %s', id)
                return None
                
            rxn.products.append(CompoundWithCoeff(coeff, compounds[id], name))
        
        return rxn
    
    @staticmethod
    def _GetCollectionAtomDiff(collection):
        """Get the net atom counts from the collection.
        
        Args:
            collection: an iterable of CompoundWithCoeff instances.
        """
        atom_diff = {}
        for compound_w_coeff in collection:
            c = compound_w_coeff.compound
            coeff = compound_w_coeff.coeff

            atom_bag = c.GetAtomBag()
            if not atom_bag:
                logging.error('Failed to fetch atom bag for %s', c.formula)
                return None
            
            for atomic_number, atom_count in atom_bag.iteritems():
                new_diff = atom_diff.get(atomic_number, 0) + coeff * atom_count
                atom_diff[atomic_number] = new_diff
        
        return atom_diff
    
    def _GetAtomDiff(self):
        """Returns the net atom counts from this reaction."""
        minus_products = [c.Minus() for c in self.products]
        return self._GetCollectionAtomDiff(self.reactants + minus_products)
    
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
        
        # Always ignore hydrogens, ala Alberty.
        atom_diff.pop('H', 0)
                        
        return max([abs(x) for x in atom_diff.values()]) < 0.01
    
    def GetBalanceWithWaterLink(self, ph=None, ionic_strength=None,
                                concentration_profile=None):
        """Returns a link to balance this reaction with water."""
        params = []
        for compound in self.reactants:
            params.append('reactantsId=%s' % compound.compound.kegg_id)
            params.append('reactantsCoeff=%d' % compound.coeff)
            if compound.name:
                params.append('reactantsName=%s' % compound.name)
        
        for compound in self.products:
            params.append('productsId=%s' % compound.compound.kegg_id)
            params.append('productsCoeff=%d' % compound.coeff)
            if compound.name:
                params.append('productsName=%s' % compound.name)
        
        if ph:
            params.append('ph=%f' % ph)
        if ionic_strength:
            params.append('ionic_strength=%f' % ionic_strength)
        if concentration_profile:
            params.append('concentration_profile=%s' % concentration_profile)
        params.append('balance_w_water=1')
    
        return '/reaction?%s' % '&'.join(params)
    
    def IsBalanced(self):
        """Checks if the collection is atom-wise balanced.
        
        Returns:
            True if the collection is atom-wise balanced.
        """
        return self._IsBalanced(self._GetAtomDiff())
    
    def _ExtraWaters(self):
        atom_diff = self._GetAtomDiff()
        if not atom_diff:
            return None
                
        # Ignore hydrogen.
        atom_diff.pop('H')
        
        # Omit oxygen for checking balancedness.
        oxy_count = atom_diff.pop('O', 0)
    
        if not self._IsBalanced(atom_diff):
            return None
        
        return oxy_count
    
    def CanBalanceWithWater(self):
        """Returns True if balanced with or without water."""
        extra_waters = self._ExtraWaters()
        if extra_waters == None:
            return False
        return True

    def TryBalanceWithWater(self):
        """Try to balance the reaction with water.
        
        Returns:
            True if the reaction is balanced already or with
            additional waters on either side.
        """ 
        extra_waters = self._ExtraWaters()
        if extra_waters == None:
            return False
        
        if extra_waters == 0:
            return True
                
        water = Compound.objects.get(kegg_id='C00001')
        w_w_coeff = CompoundWithCoeff(compound=water, coeff=abs(extra_waters),
                                      name='Water')
        if extra_waters > 0:
            self.products.append(w_w_coeff)
        else:
            self.reactants.append(w_w_coeff)
        
        return True

    @staticmethod
    def GetTotalFormationEnergy(collection,
                                pH=constants.DEFAULT_PH,
                                ionic_strength=constants.DEFAULT_IONIC_STRENGTH):
        """Compute an estimate for a collection of compounds + coefficients.
        
        You can compute the DeltaG of a reaction using negative coefficients for
        products.
        
        Args:
            collection: an iterable of CompoundWithCoeff objects.
        """        
        sum = 0
        for compound_w_coeff in collection:
            c = compound_w_coeff.compound
            coeff = compound_w_coeff.coeff
            
            est = c.DeltaG(pH, ionic_strength)
            if not est:
                logging.info('No estimate for compound %s', id)
                return None
            
            sum += coeff * est
        
        return sum
    
    def _GetConcentrationCorrection(self, concentration_profile):
        """Get the concentration term in DeltaG' for these concentrations.
        
        Args:
            concentration_profile: a ConcentrationProfile object.
        
        Returns:
            The correction or None on error.
        """
        reactant_term, product_term = 0, 0
        cp = concentration_profile
        
        try:
            mult_log = lambda c: c.coeff * pylab.log(cp.Concentration(c.compound.kegg_id))
            reactant_terms = [mult_log(c) for c in self.reactants]
            product_terms = [mult_log(c) for c in self.products]
            reactant_term = sum(reactant_terms)
            product_term = sum(product_terms)
        except KeyError, e:
            logging.error(e)
            return None
        
        _r = constants.R
        _t = constants.DEFAULT_TEMP
        return _r * _t * (product_term - reactant_term)
    
    def DeltaG(self,
               pH=constants.DEFAULT_PH,
               ionic_strength=constants.DEFAULT_IONIC_STRENGTH,
               concentration_profile=None):
        """Compute the DeltaG for a reaction.
        
        Args:
            pH: the PH to estimate at.
            ionic_strength: the ionic strength to estimate at.
            temp: the temperature to estimate at.
            concentrations: a dictionary mapping kegg IDs to concentrations.
        
        Returns:
            The DeltaG for this reaction, or None if data was missing.
        """
        reactants_sum = self.GetTotalFormationEnergy(
            self.reactants, pH, ionic_strength)
        products_sum = self.GetTotalFormationEnergy(
            self.products, pH, ionic_strength)
        if not products_sum:
            logging.warning('Failed to get products formation energy.')
            return None
        if not reactants_sum:
            logging.warning('Failed to get reactants formation energy.')
            return None
        
        dg_zero = products_sum - reactants_sum
        if not concentration_profile:
            return dg_zero        
        return dg_zero + self._GetConcentrationCorrection(concentration_profile)
    
    def NoDeltaGExplanation(self):
        """Get an explanation for why there's no delta G value.
        
        Return:
            The explanation or None.
        """
        for compound in self.reactants + self.products:
            if compound.compound.no_dg_explanation:
                name = compound.compound.common_names.all()[0].name
                return '%s %s' % (name,
                                  compound.compound.no_dg_explanation.lower())
        return None
        
        
    
