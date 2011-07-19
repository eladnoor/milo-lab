#!/usr/bin/python

import hashlib
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
    
    # A citation name for the source. May be null.
    citation = models.CharField(max_length=4096, null=True)
    
    # The pubmed ID of the source if available.
    pubmed_id = models.CharField(max_length=128, null=True)
    
    # A link explaining the source.
    link = models.URLField(null=True)
    
    def __unicode__(self):
        return self.name
    
    @staticmethod
    def _GetOrCreate(name):
        """Gets or creates a ValueSource object with the given name from the DB.
        
        TODO(flamholz): supply a link here?
        
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
    def Thauer():
        return ValueSource._GetOrCreate('Thauer et al.')
    
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
    formation_energy_source = models.ForeignKey(ValueSource, null=True)
    
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
    
class SpeciesGroup(models.Model):
    # The ID of the compound in KEGG.
    kegg_id = models.CharField(max_length=10)
    
    # The species in this group.
    species = models.ManyToManyField(Specie)
    
    # The priority of this group.
    priority = models.IntegerField()
    
    # The source of these values.
    formation_energy_source = models.ForeignKey(ValueSource)
    
    def __init__(self, *args, **kwargs):        
        super(SpeciesGroup, self).__init__(*args, **kwargs)
        self._all_species = None
            
    def GetSpecies(self):
        """Gets the list of SpeciesFormationEnergies, potentially caching."""
        if not self._all_species:
            self._all_species = self.species.all()
        return self._all_species
    all_species = property(GetSpecies)

    def StashTransformedSpeciesEnergies(self, ph, pmg, ionic_strength):
        """Stash the transformed species formation energy in each one."""
        for species in self.all_species:
            species.transformed_energy = species.Transform(
                pH=ph, pMg=pmg, ionic_strength=ionic_strength)

    def DeltaG(self, pH=constants.DEFAULT_PH,
               pMg=constants.DEFAULT_PMG,
               ionic_strength=constants.DEFAULT_IONIC_STRENGTH):
        """Get a deltaG estimate for this group of species.
        
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


class Compound(models.Model):
    """A single compound."""
    # The ID of the compound in KEGG.
    kegg_id = models.CharField(max_length=10, null=True)
    
    # InChI representation of the compound.
    inchi = models.CharField(max_length=2048, null=True)
        
    # A list of common names of the compound, used for searching.
    common_names = models.ManyToManyField(CommonName)
    
    # If present, this name should always be used.
    preferred_name = models.CharField(max_length=250, null=True)
    
    # A remark about this compound.
    note = models.TextField(null=True)
    
    # A link for detailed remarks about this compound. 
    details_link = models.URLField(null=True)

    # The chemical formula.
    formula = models.CharField(max_length=500, null=True)
    
    # The molecular mass.
    mass = models.FloatField(null=True)  # In Daltons.
    
    # The number of electrons.
    num_electrons = models.IntegerField(null=True)
    
    # The various estimates of Delta G for this compound.
    species_groups = models.ManyToManyField(SpeciesGroup)
    
    # Replace this compound with another one.
    replace_with = models.ForeignKey('self', null=True)
    
    # An explanation for when no DeltaG estimate is available.
    no_dg_explanation = models.CharField(max_length=2048,
                                         blank=True,
                                         null=True)

    def __init__(self, *args, **kwargs):        
        super(Compound, self).__init__(*args, **kwargs)
        self._species_group_to_use = None
        self._priority = None
    
    def GetSpeciesGroupPriorities(self):
        """Get the priorities of species groups available."""
        return [sg.priority for sg in self.all_species_groups]
    
    def SetSpeciesGroupPriority(self, priority):
        """Set the priority of the species group to use."""
        for sg in self.species_groups.all():
            if sg.priority == priority:
                self._species_group_to_use = sg
                break
    
    def SetHighestPriority(self):
        """Set the priority to the highest one."""
        ps = self.GetSpeciesGroupPriorities()
        if not ps:
            return
        
        self.SetSpeciesGroupPriority(min(ps))
    
    def HasData(self):
        """Has enough data to display."""
        return self.mass and self.formula
    
    def FirstName(self):
        """Return the 'first' name of this compound.
        
        If a 'preferred_name' is set, returns that. Otherwise, returns
        the first name in the list of common names. Presumes that the
        list of names is in some order.
        """
        if self.preferred_name:
            return self.preferred_name
        
        names = list(self.common_names.all())
        return names[0].name
    
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
        return self._species_group.DeltaG(pH, pMg, ionic_strength)

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
        """Returns the link to the stand-alone page for this compound."""
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
    
    def GetSpeciesGroupToUse(self):
        """Gets the SpeciesGroup to use, potentially caching."""
        if self._species_group_to_use:
            return self._species_group_to_use
        
        self.SetHighestPriority()
        return self._species_group_to_use
    
    def GetSpecies(self):
        """Gets the list of SpeciesFormationEnergies, potentially caching."""
        if self._species_group_to_use:
            return self._species_group_to_use.all_species
        
        # TODO(flamholz): Should we return something here?
        return None

    def GetSpeciesGroups(self):
        """Gets the list of SpeciesGroups."""
        return self.species_groups.all()

    def _GetDGSource(self):
        """Returns the source of the dG data."""
        if not self._species_group:
            return None
        
        return self._species_group.formation_energy_source        
    
    def _GetSubstrateEnzymes(self):
        return self.substrate_for_enzymes.all()
    
    def _GetProductEnzymes(self):
        return self.product_of_enzymes.all()
    
    def _GetCofactorEnzymes(self):
        return self.cofactor_of_enzymes.all()
    
    _species_group = property(GetSpeciesGroupToUse)
    first_name = property(FirstName)
    html_formula = property(GetHtmlFormattedFormula)
    link = property(GetLink)
    kegg_link = property(GetKeggLink)
    small_image_url = property(GetSmallImageUrl)
    all_common_names = property(lambda self: self.common_names.all())
    all_species = property(GetSpecies)
    all_species_groups = property(GetSpeciesGroups)
    substrate_of = property(_GetSubstrateEnzymes)
    product_of = property(_GetProductEnzymes)
    cofactor_of = property(_GetCofactorEnzymes)
    standard_formation_energy = property(DeltaG)
    dg_source = property(_GetDGSource)
    
    def StashTransformedSpeciesEnergies(self, ph, pmg, ionic_strength):
        """Stash the transformed species formation energy in each one."""
        for sg in self.all_species_groups:
            sg.StashTransformedSpeciesEnergies(ph, pmg, ionic_strength)
    
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
    

class Reactant(models.Model):
    """A compound and its coefficient."""
    # The compound.
    compound = models.ForeignKey(Compound)
    
    # The coeff.
    coeff = models.IntegerField(default=1)

    @staticmethod
    def GetOrCreate(kegg_id, coeff):
        """Attempts to fetch the CommonName object, or creates it if not present.
        
        Args:
            name: the name to use.
        
        Returns:
            A CommonName object.
        """
        try:
            r = Reactant.objects.get(compound__kegg_id=kegg_id,
                                     coeff=coeff)
            return r
        except Exception:
            c = Compound.objects.get(kegg_id=kegg_id)
            r = Reactant(compound=c, coeff=coeff)
            r.save()
            return r
    

class StoredReaction(models.Model):
    """A reaction stored in the database."""
    # The ID of this reaction in KEGG.
    kegg_id = models.CharField(max_length=10, null=True)
    
    # The list of reactants.
    reactants = models.ManyToManyField(Reactant, related_name='reactant_in')
    
    # The list of products.
    products = models.ManyToManyField(Reactant, related_name='product_in')
    
    # TODO(flamholz): add some sort of hash identifier.
    hash = models.CharField(max_length=2048)
    
    @staticmethod
    def _SideString(side):
        """Returns a string representation for a single side of a reaction.
        
        Args:
            side: the list of CompoundWithCoeff objects representing the side.
        """
        l = []
        for r in side:
            if r.coeff == 1:
                l.append(r.compound.FirstName())
            else:
                l.append('%d %s' % (r.coeff,
                                    r.compound.FirstName()))
        return ' + '.join(l)

    def ReactionString(self):
        """Get the string representation of this reaction."""
        return '%s <=> %s' % (self._SideString(self.reactants.all()),
                              self._SideString(self.products.all()))
    
    @staticmethod
    def HashableReactionString(reactants, products):
        """Return a hashable string for a biochemical reaction.
        
        The string fully identifies the biochemical reaction up to directionality.
        If it is equal to another reaction's string, then they have identical
        stoichiometry up to their directionality.
        
        Args:
            reactants: the reactants; a list of Reactants or like objects.
            products: the products; a list of Reactants or like objects.
        """
        sort_key = lambda r: r.compound.kegg_id
        make_str = lambda r: '%d %s' % (r.coeff, r.compound.kegg_id)
        is_not_hydrogen = lambda r: r.compound.kegg_id != 'C00080'
        
        reactants_strs = map(make_str,
                             sorted(filter(is_not_hydrogen, reactants),
                                    key=sort_key))
        rside_str = ' + '.join(reactants_strs)
        rside_hash = str(hash(rside_str))
        
        products_strs = map(make_str,
                            sorted(filter(is_not_hydrogen, products),
                                   key=sort_key))
        pside_str = ' + '.join(products_strs)
        pside_hash = str(hash(pside_str))
        
        sides = ['%s%s' % (rside_hash, rside_str),
                 '%s%s' % (pside_hash, pside_str)]
        sides.sort()
        return '%s <=> %s' % (sides[0], sides[1])
    
    def GetHashableReactionString(self):
        """Get a hashable string identifying this chemical reaction."""
        return self.HashableReactionString(self.reactants.all(),
                                           self.products.all())
    
    @staticmethod
    def HashReaction(reactants, products):
        md5 = hashlib.md5()
        md5.update(StoredReaction.HashableReactionString(reactants, products))
        return md5.hexdigest()
    
    def GetHash(self):
        """Returns a string hash of this reaction for easy identification."""
        return self.HashReaction(self.reactants.all(),
                                 self.products.all())
    
    def __hash__(self):
        """Makes stored reactions hashable."""
        return hash(self.GetHash())
    
    def __str__(self):
        """String representation."""
        return self.ReactionString()
    
    def Link(self):
        """Returns a link to this reaction's page."""
        return '/reaction?reactionId=%s' % self.kegg_id
    
    link = property(Link)
    reaction_string = property(ReactionString)


class Enzyme(models.Model):
    """A single enzyme."""
    # EC class enzyme.
    ec = models.CharField(max_length=10)
    
    # A list of common names of the compound, used for searching.
    common_names = models.ManyToManyField(CommonName)
    
    # List of reactions this enzyme catalyzes.
    reactions = models.ManyToManyField(StoredReaction)
    
    # Compounds that this reaction
    substrates = models.ManyToManyField(Compound, related_name='substrate_for_enzymes')
    products = models.ManyToManyField(Compound, related_name='product_of_enzymes')
    cofactors = models.ManyToManyField(Compound, related_name='cofactor_of_enzymes')
    
    # TODO(flamholz): add more fields.
    
    def HasData(self):
        """Checks if it has enough data to display."""
        return self.ec and self.reactions.all()

    def Link(self):
        """Returns a link to this reactions page."""
        return '/enzyme?ec=%s' % self.ec
    
    def KeggLink(self):
        """Returns a link to the KEGG page for this enzyme."""
        return 'http://kegg.jp/dbget-bin/www_bget?ec:%s' % self.ec
    
    def BrendaLink(self):
        """Returns a link to the BRENDA page for this enzyme."""
        return 'http://www.brenda-enzymes.org/php/result_flat.php4?ecno=%s' % self.ec

    def AllReactions(self):
        """Returns all the reactions."""
        return self.reactions.all()

    def FirstName(self):
        """The first name in the list of names."""
        return self.all_common_names[0]

    def __unicode__(self):
        """Return a single string identifier of this enzyme."""
        return unicode(self.FirstName())

    all_common_names = property(lambda self: self.common_names.all())
    all_cofactors = property(lambda self: self.cofactors.all())
    first_name = property(FirstName)
    all_reactions = property(AllReactions)
    kegg_link = property(KeggLink)
    brenda_link = property(BrendaLink)
    link = property(Link)
