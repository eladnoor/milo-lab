import logging
import numpy
import re

from django.db import models
from gibbs import constants
from util import inchi
        

class CommonName(models.Model):
    """A common name of a compound."""
    name = models.CharField(max_length=500)
    enabled = models.BooleanField()
    
    @staticmethod
    def GetOrCreate(name):
        """Attempts to fetch the CommonName object, or creates it if not present.
        
        Args:
            name: the name to use.
        
        Returns:
            A CommonName object.
        """
        try:
            n = CommonName.objects.get(name=name)
            return n
        except Exception:
            n = CommonName(name=name)
            n.save()
            return n
    
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
    
    @staticmethod
    def _GetOrCreate(name):
        """Gets or creates a ValueSource object with the given name from the DB.
        
        Args:
            name: the name field.
        
        Returns:
            A ValueSource object.
        """
        vs = None
        try:
            vs = ValueSource.objects.get(name=name)
        except Exception:
            vs = ValueSource(name=name)
            vs.save()
        return vs
    
    @staticmethod
    def Alberty():
        return ValueSource._GetOrCreate('Alberty et al.')
    
    @staticmethod
    def GroupContribution():
        return ValueSource._GetOrCreate('Estimated using group contribution')

class Specie(models.Model):
    """A single specie of a compound."""
    # The ID of the compound in KEGG.
    kegg_id = models.CharField(max_length=10)
    
    # The number of hydrogens in the species.
    number_of_hydrogens = models.IntegerField()
    
    # The number of Mg2+ ions bound to this species.
    number_of_mgs = models.IntegerField(default=0)
    
    # The net charge (eV).
    net_charge = models.IntegerField()
    
    # The standard formation energy in kJ/mol.
    formation_energy = models.FloatField()
    
    # The source of this value.
    formation_energy_source = models.ForeignKey(ValueSource)
    
    def Transform(self,
                  pH=constants.DEFAULT_PH,
                  pMg=constants.DEFAULT_PMG,
                  ionic_strength=constants.DEFAULT_IONIC_STRENGTH):
        """Transform this individual estimate to difference conditions."""
        # Short names are nice!
        _i_s = ionic_strength
        _r = constants.R
        _t = constants.DEFAULT_TEMP
        _dgf_mg = constants.MG_FORMATION_ENERGY
        _n_h = self.number_of_hydrogens
        _n_mg = self.number_of_mgs
        _n_c = self.net_charge
        _dg = self.formation_energy
        
        chem_potential = _n_h * _r * _t * numpy.log(10) * pH
        ionic_potential = (2.91482 * (_n_c ** 2 - _n_h) * numpy.sqrt(_i_s) /
                           (1 + 1.6 * numpy.sqrt(_i_s)))
        mg_potential = _n_mg * (_r * _t * numpy.log(10) * pMg - _dgf_mg)
        return _dg + mg_potential + chem_potential - ionic_potential
    
    def __unicode__(self):
        return self.kegg_id


class Compound(models.Model):
    """A single compound."""
    # The ID of the compound in KEGG.
    kegg_id = models.CharField(max_length=10, null=True)
    
    # InChI representation of the compound.
    inchi = models.CharField(max_length=2048, null=True)
    
    # InChI string with chirality fields removed.
    achiral_inchi = models.CharField(max_length=2048, null=True)
    
    # A list of common names of the compound, used for searching.
    common_names = models.ManyToManyField(CommonName)
    
    # The chemical formula.
    formula = models.CharField(max_length=500, null=True)
    
    # The molecular mass.
    mass = models.FloatField(null=True)  # In Daltons.
    
    # The number of electrons.
    num_electrons = models.IntegerField(null=True)
    
    # Estimates of Delta G for this compound.
    species = models.ManyToManyField(Specie)
    
    # An explanation for when no DeltaG estimate is available.
    no_dg_explanation = models.CharField(max_length=2048,
                                         blank=True,
                                         null=True)

    def __init__(self, *args, **kwargs):        
        super(Compound, self).__init__(*args, **kwargs)
        self._all_species = None
    
    def save(self):
        """Custom save-time behavior."""
        if self.inchi:
            self.achiral_inchi = inchi.AchiralInchi(self.inchi)
        
        super(Compound, self).save()
    
    def clean(self):
        """Custom save-time behavior."""
        if self.inchi:
            self.achiral_inchi = inchi.AchiralInchi(self.inchi)
    
    def FirstName(self):
        return self.common_names.all()[0]
    
    def ShortestName(self):
        shortest_len = 10000
        shortest_name = None
        for name in self.common_names.all():
            if len(name.name) < shortest_len:
                shortest_name = name.name
                shortest_len = len(name.name)
                
        return shortest_name
    
    def DeltaG(self, pH=constants.DEFAULT_PH,
               pMg=constants.DEFAULT_PMG,
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
        transform = lambda x: x.Transform(pH=pH, pMg=pMg, ionic_strength=_i_s)
        scaled_transforms = [(-transform(s) / (_r * _t))
                             for s in self.all_species]
        
        # Numerical issues: taking a sum of exp(v) for |v| quite large.
        # Use the fact that we take a log later to offset all values by a 
        # constant (the minimum value).
        offset = min(scaled_transforms)
        scaled_offset_transforms = [(st - offset) for st in scaled_transforms]
        sum_exp = sum(numpy.exp(scaled_offset_transforms))
        return - _r * _t * (offset + numpy.log(sum_exp))

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
    
    def _GetDGSource(self):
        """Returns the source of the dG data."""
        if not self.all_species:
            return None
        
        return self.all_species[0].formation_energy_source
        
    html_formula = property(GetHtmlFormattedFormula)
    link = property(GetLink)
    kegg_link = property(GetKeggLink)
    small_image_url = property(GetSmallImageUrl)
    all_common_names = property(lambda self: self.common_names.all())
    all_species = property(GetSpecies)
    standard_formation_energy = property(DeltaG)
    dg_source = property(_GetDGSource)
    
    def StashTransformedSpeciesEnergies(self, ph, ionic_strength):
        """Stash the transformed species formation energy in each one."""
        for species in self.all_species:
            species.transformed_energy = species.Transform(
                pH=ph, ionic_strength=ionic_strength)
    
    def __unicode__(self):
        """Return a single string identifier of this Compound."""
        names = self.all_common_names
        if names:
            return unicode(names[0])
        return unicode(self.formula)
    
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

