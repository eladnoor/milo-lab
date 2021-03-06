# -*- coding: utf-8 -*-

import logging
import numpy
import urllib

from gibbs import concentration_profile
from gibbs import constants
from gibbs import models

# gases by order in KEGG: O2, NH3, CO, H2, N2   
VOLATILE_COMPOUNDS_KEGG_IDS = ['C00007','C00014','C00237','C00282','C00697']

class ReactantFormulaMissingError(Exception):
    
    def __init__(self, c):
        self.compound = c
        
    def __str__(self):
        return "Cannot test reaction balancing because the reactant %s does not have a chemical formula" % self.compound.kegg_id

class CompoundWithCoeff(object):
    """A compound with a stoichiometric coefficient."""
    
    def __init__(self, coeff, compound, name=None, concentration=1.0):
        """Construct a CompoundWithCoeff object.
        
        Args:
            coeff: the coefficient.
            compound: the Compound object.
            name: a string of the compound name.
            concentration: the concentration (molar).
        """
        self.compound = compound
        self.coeff = coeff
        self._name = name
        self.concentration = concentration
        self.transformed_energy = None
        
    @staticmethod
    def FromReactant(reactant):
        return CompoundWithCoeff(reactant.coeff, reactant.compound,
                                 name=reactant.compound.FirstName())
    
    @staticmethod
    def FromId(coeff, kegg_id, name=None, concentration=1.0):
        compound = models.Compound.objects.get(kegg_id=kegg_id)
        my_name = name or compound.FirstName()
        return CompoundWithCoeff(coeff, compound, name=my_name,
                                 concentration=concentration)
        
    def Minus(self):
        """Returns a new CompoundWithCoeff with coeff = -self.coeff."""
        return CompoundWithCoeff(-self.coeff, self.compound, self.name)
    
    def __str__(self):
        name = self.name or str(self.compound)
        return '%d %s' % (self.coeff, name)
    
    def ToJson(self, include_species=True):
        d = {'coeff': self.coeff,
             'KEGG_ID': self.compound.kegg_id,
             'concentration': self.concentration,
             'name': str(self.compound.first_name),
             'source_used': None}
        
        if self.compound._species_group is not None:
            d['source_used'] = str(self.compound._species_group.formation_energy_source)
        
        if include_species:
            d["species"] = self.compound.SpeciesJson()
            
        return d
    
    def GetName(self):
        """Gives a string name for this compound."""
        if self.compound.preferred_name:
            return self.compound.preferred_name
        if self._name:
            return self._name
        return str(self.compound.FirstName())
    
    def _MicromolarConcentration(self):
        return self.concentration * 1e6
    
    def _MicromolarConcentrationString(self):
        conc = self.concentration * 1e6
        if conc > 1000:
            return '%.2e' % conc
        return '%.2f' % conc
    
    def _HumanConcentrationStringWithUnits(self):
        if self.concentration > 1e-2:
            return '%.2g M' % self.concentration
        
        if self.concentration > 1e-4:
            return '%.2g mM' % (self.concentration * 1e3)
        
        return '%2g μM' % (self.concentration * 1e6)
    
    def __eq__(self, other):
        """Check equality with another CompoundWithCoeff.
        
        Args:
            other: a second CompoundWithCoeff (or like object).
        """
        if self.coeff != other.coeff:
            return False
        
        if self.compound.kegg_id != other.compound.kegg_id:
            return False
        
        return True
    
    name = property(GetName)
    micromolar_concentration = property(_MicromolarConcentration)
    micromolar_concentration_string = property(_MicromolarConcentrationString)
    human_concentration_w_units = property(_HumanConcentrationStringWithUnits)
    

class Reaction(object):
    """A reaction."""
    
    def __init__(self, substrates=None, products=None,
                 pH=constants.DEFAULT_PH,
                 pMg=constants.DEFAULT_PMG,
                 ionic_strength=constants.DEFAULT_IONIC_STRENGTH):
        """Construction.
        
        Args:
            substrates: a list of CompoundWithCoeff objects.
            products: a list of CompoundWithCoeff objects.
        """
        self.substrates = self._FilterZeroes(substrates or [])
        self.products = self._FilterZeroes(products or [])
        self.ph = pH
        self.pmg = pMg
        self.i_s = ionic_strength
        self._concentration_profile = None
        self._stored_reaction = None
        self._all_stored_reactions = None
        self._catalyzing_enzymes = None
        self._SetCompoundPriorities()
    
    def _SetCompoundPriorities(self):
        """Returns a set of (int, SpeciesGroup) tuples for the reaction."""
        # Hack - ignore H+ because it has the wrong priority for the Alberty data.
        compounds = self._FilterHydrogen(self.substrates + self.products)
        compounds = [c.compound for c in compounds]
        
        sentinel = 1<<10;
        get_min_priority = lambda c: min(c.GetSpeciesGroupPriorities() + [sentinel])
        min_priorities = map(get_min_priority, compounds)
        priority_to_use = max(min_priorities)
        
        # Someone is missing data!
        if priority_to_use == sentinel:
            return
        
        for c in compounds:
            c.SetSpeciesGroupPriority(priority_to_use)
         
    def SwapSides(self):
        """Swap the sides of this reaction."""
        tmp = self.substrates
        self.substrates = self.products
        self.products = tmp
        
    def GetConcentrationProfile(self):
        """Get the concentration profile of this reaction."""
        return self._concentration_profile
    
    def ApplyConcentrationProfile(self, concentration_profile):
        """Apply this concentration profile to this reaction.
        
        Args:
            concentration_profile: a ConcentrationProfile object.
        """
        self._concentration_profile = concentration_profile
        for c in self.substrates + self.products:
            c.concentration = self._concentration_profile.Concentration(c.compound.kegg_id)
    concentration_profile = property(GetConcentrationProfile,
                                     ApplyConcentrationProfile)    
    
    def __str__(self):
        """Simple text reaction representation."""
        rlist = map(str, self.substrates)
        plist = map(str, self.products)
        return '%s <=> %s' % (' + '.join(rlist), ' + '.join(plist))
    
    def SameChemicalReaction(self, stored_reaction):
        """Checks that the two chemical reactions are the same."""
        my_string = models.StoredReaction.HashableReactionString(self.substrates,
                                                                 self.products)
        their_string = stored_reaction.GetHashableReactionString()
        return my_string == their_string
    
    def _GetHash(self):
        """Get the hash for this reaction."""
        if self._stored_reaction:
            return self._stored_reaction.GetHash()
        
        return models.StoredReaction.HashReaction(self.substrates,
                                                  self.products)
    
    def _GetAllStoredReactions(self):
        """Find all stored reactions matching this compound."""
        if not self._all_stored_reactions:
            self._all_stored_reactions = []
            hash = self._GetHash()
            matching_reactions = models.StoredReaction.objects.select_related(
                ).filter(hash=hash)
            self._all_stored_reactions = filter(self.SameChemicalReaction,
                                                matching_reactions)
        
        return self._all_stored_reactions
    all_stored_reactions = property(_GetAllStoredReactions)
    
    def GetStoredReaction(self):
        """Get the stored reaction if any."""
        return self._stored_reaction
    
    def SetStoredReaction(self, stored_reaction):
        """Setter for stored reaction."""
        self._stored_reaction = stored_reaction
    stored_reaction = property(GetStoredReaction,
                               SetStoredReaction)
    
    def GetCatalyzingEnzymes(self):
        """Get all the enzymes catalyzing this reaction."""        
        if not self._catalyzing_enzymes:
            self._catalyzing_enzymes = []
            for stored_reaction in self.all_stored_reactions:
                enzymes = stored_reaction.enzyme_set.all()
                self._catalyzing_enzymes.extend(enzymes)
        
        return self._catalyzing_enzymes
    catalyzing_enzymes = property(GetCatalyzingEnzymes)
    
    def StandardConcentrations(self):
        """Returns True if using standard concentrations."""
        if not self._concentration_profile:
            return True
        return (self._concentration_profile and
                self._concentration_profile.IsStandard())
    
    def ToJson(self):
        """Return this reaction as a JSON-compatible object."""
        pdicts = [c.ToJson(include_species=False) for c in self.products]
        rdicts = [r.ToJson(include_species=False) for r in self.substrates]
        enzdicts = [e.ToJson() for e in self.catalyzing_enzymes]
        d = {'reaction_string': str(self), 
             'substrates': rdicts,
             'products': pdicts,
             'enzymes': enzdicts,
             'chemically_balanced': self.is_balanced,
             'redox_balanced': self.is_electron_balanced,
             'dgzero': None,
             'dgzero_tag': None,
             'keq_tag': None,
             'KEGG_ID': None}
        
        if self.dg0 is not None:
            d['dgzero'] = round(self.dg0, 1)
        if self.dg0_tag is not None:
            d['dgzero_tag'] = {
                'value': round(self.dg0_tag, 1),
                'pH': self.ph,
                'ionic_strength': self.i_s}
        if self.k_eq_tag is not None:
            d['keq_tag'] = {
                'value': self.k_eq_tag,
                'pH': self.ph,
                'ionic_strength': self.i_s}
        
        if self.stored_reaction:
            d['KEGG_ID'] = self.stored_reaction.kegg_id
        
        return d
    
    @staticmethod
    def FromForm(form):
        """Build a reaction object from a ReactionForm.
        
        Args:
            form: a ReactionForm object.
        
        Returns:
            A Reaction object or None if there's an error.
        """
        if form.cleaned_reactionId:
            stored_reaction = models.StoredReaction.objects.get(
                kegg_id=form.cleaned_reactionId)
            return Reaction.FromStoredReaction(stored_reaction)
        
        i_s = form.cleaned_ionic_strength
        ph = form.cleaned_ph
        pmg = form.cleaned_pmg
        
        clean_substrates = form.cleaned_substrateIds
        clean_products = form.cleaned_productIds
        all_ids = clean_substrates + clean_products
        
        # Fetch custom concentrations if any.
        substrate_concentrations = form.cleaned_substrateConcentrations
        product_concentrations = form.cleaned_productConcentrations
        all_concentrations = substrate_concentrations + product_concentrations
        
        # Build the appropriate concentration profile.
        cprofile_name = form.cleaned_concentration_profile
        cprofile = concentration_profile.GetProfile(
            cprofile_name, all_ids, all_concentrations)
        
        substrate_names = form.cleaned_substrateNames
        product_names = form.cleaned_productNames
        
        # Return the built reaction object.
        zipped_substrates = zip(form.cleaned_substrateCoeffs, clean_substrates, substrate_names)
        zipped_products = zip(form.cleaned_productCoeffs, clean_products, product_names)
        return Reaction.FromIds(zipped_substrates, zipped_products,
                                concentration_profile=cprofile,
                                pH=ph, pMg=pmg,
                                ionic_strength=i_s)
    
    @staticmethod
    def FromStoredReaction(stored_reaction,
                           pH=constants.DEFAULT_PH,
                           pMg=constants.DEFAULT_PMG,
                           ionic_strength=constants.DEFAULT_IONIC_STRENGTH):
        """Build a reaction object from a stored reaction.
        
        Args:
            stored_reaction: models.StoredReaction object.
            
        Returns:
            A Reaction object.
        """
        substrates = [CompoundWithCoeff.FromReactant(r)
                     for r in stored_reaction.substrates.all()]
        products  = [CompoundWithCoeff.FromReactant(r)
                     for r in stored_reaction.products.all()]
        rxn = Reaction(substrates, products, pH=pH, pMg=pMg,
                       ionic_strength=ionic_strength)
        rxn.SetStoredReaction(stored_reaction)
        return rxn
    
    @staticmethod
    def FromIds(substrates, products, concentration_profile=None,
                pH=constants.DEFAULT_PH,
                pMg=constants.DEFAULT_PMG,
                ionic_strength=constants.DEFAULT_IONIC_STRENGTH):
        """Build a reaction object from lists of IDs.
        
        Args:
            substrates: an iterable of (coeff, kegg_id, name) of substrates.
            products: an iterable of (coeff, kegg_id, name) of products.
            concentration_profile: a ConcentrationProfile object.
            
        Returns:
            A properly set-up Reaction object or None if there's an error.
        """
        r_ids = [id for unused_coeff, id, unused_name in substrates]
        p_ids = [id for unused_coeff, id, unused_name in products]
        compounds = models.Compound.GetCompoundsByKeggId(r_ids + p_ids)
        
        # Get products and substrates.
        rs, ps = [], []
        for coeff, id, name in substrates:
            if id not in compounds:
                logging.error('Unknown substrate %s', id)
                return None
            
            rs.append(CompoundWithCoeff(coeff, compounds[id], name))
        
        for coeff, id, name in products:
            if id not in compounds:
                logging.error('Unknown product %s', id)
                return None
                
            ps.append(CompoundWithCoeff(coeff, compounds[id], name))
        
        rxn = Reaction(rs, ps, pH=pH, pMg=pMg,
                       ionic_strength=ionic_strength)
        if concentration_profile:
            rxn.ApplyConcentrationProfile(concentration_profile)
        
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
                logging.warning('Failed to fetch atom bag for %s', c.formula)
                raise ReactantFormulaMissingError(c)
            
            for atomic_number, atom_count in atom_bag.iteritems():
                new_diff = atom_diff.get(atomic_number, 0) + coeff * atom_count
                atom_diff[atomic_number] = new_diff
        
        return atom_diff
    
    @staticmethod
    def _GetCollectionElectronDiff(collection):
        """Get the net electron count from the collection.
        
        Args:
            collection: an iterable of CompoundWithCoeff instances.
        """
        electron_diff = 0
        for compound_w_coeff in collection:
            c = compound_w_coeff.compound
            coeff = compound_w_coeff.coeff
            
            electrons = c.num_electrons
            if electrons == None:
                logging.warning('Compound %s has unknown electron count', c.kegg_id)
                return 0
            
            electron_diff += coeff * electrons
        return electron_diff
    
    def _GetAtomDiff(self):
        """Returns the net atom counts from this reaction."""
        minus_products = [c.Minus() for c in self.products]
        return self._GetCollectionAtomDiff(self.substrates + minus_products)
    
    def _GetElectronDiff(self):
        """Returns the net electron count from this reaction."""
        minus_products = [c.Minus() for c in self.products]
        return self._GetCollectionElectronDiff(self.substrates + minus_products)
    
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
    
    def _GetUrlParams(self, query=None):
        """Get the URL params for this reaction."""
        params = []
        for compound in self.substrates:
            params.append('substratesId=%s' % compound.compound.kegg_id)
            params.append('substratesCoeff=%d' % compound.coeff)
            if compound.name:
                params.append('substratesName=%s' % compound.name)
        
        for compound in self.products:
            params.append('productsId=%s' % compound.compound.kegg_id)
            params.append('productsCoeff=%d' % compound.coeff)
            if compound.name:
                params.append('productsName=%s' % compound.name)
        
        if self.ph:
            params.append('ph=%f' % self.ph)
        if self.pmg:
            params.append('pmg=%f' % self.pmg)
        if self.i_s:
            params.append('ionic_strength=%f' % self.i_s)
        if self._concentration_profile:
            params.append('concentration_profile=%s' % self._concentration_profile)
            if self._concentration_profile.IsCustom():
                l = ['substratesConcentration=%s' % self._concentration_profile.MicroMolarConcentration(c.compound.kegg_id)
                     for c in self.substrates]
                params.extend(l)
                
                l = ['productsConcentration=%s' % self._concentration_profile.MicroMolarConcentration(c.compound.kegg_id)
                     for c in self.products]
                params.extend(l)
                
        if query:
            tmp_query = query.replace(u'→', '=>')
            params.append('query=%s' % urllib.quote(tmp_query))
            
        return params
    
    def GetBalanceWithWaterLink(self, query=None):
        """Returns a link to balance this reaction with water."""
        params = self._GetUrlParams(query)
        params.append('balance_w_water=1')
    
        return '/reaction?%s' % '&'.join(params)

    def GetBalanceElectronsLink(self, query=None):
        """Returns a link to balance this reaction with NAD:NADH."""
        params = self._GetUrlParams(query)
        params.append('balance_electrons=1')
    
        return '/reaction?%s' % '&'.join(params)

    def GetHalfReactionLink(self, query=None):
        """Returns a link to balance this reaction with water."""
        params = self._GetUrlParams(query)
        return '/half_reaction?%s' % '&'.join(params)

    def GetReplaceCO2Link(self, query=None):
        """Returns a link to balance this reaction with water."""
        params = self._GetUrlParams(query)
        params.append('replace_co2=1')
    
        return '/reaction?%s' % '&'.join(params)
    
    def GetPhGraphLink(self):
        params = self._GetUrlParams()
        params.append('vary_ph=1')
        return '/graph_reaction?%s' % '&'.join(params)

    def GetPMgGraphLink(self):
        params = self._GetUrlParams()
        params.append('vary_pmg=1')
        return '/graph_reaction?%s' % '&'.join(params)
    
    def GetIonicStrengthGraphLink(self):
        params = self._GetUrlParams()
        params.append('vary_is=1')
        return '/graph_reaction?%s' % '&'.join(params)
    
    @staticmethod
    def _GetReactionSideString(side):
        """Write a reaction side as a string."""
        sdata = []
        for c_w_coeff in side:
            if c_w_coeff.coeff == 1:
                sdata.append(c_w_coeff.GetName())
            else:
                sdata.append('%d %s' % (c_w_coeff.coeff,
                                        c_w_coeff.GetName()))
        return ' + '.join(sdata)
    
    def GetQueryString(self):
        """Get a query string for this reaction."""
        rq = self._GetReactionSideString(self.substrates)
        pq = self._GetReactionSideString(self.products)
        return '%s <=> %s' % (rq, pq)
    
    def ContainsCO2(self):
        co2_id = 'C00011'
        if self._FindCompoundIndex(self.substrates, co2_id) is not None:
            return True
        return self._FindCompoundIndex(self.products, co2_id) is not None

    def ContainsVolatile(self):
        """Checks if at least one of the reactants is volatile
        
        Returns:
            True if there is a volatile reactant
        """
        for v_id in VOLATILE_COMPOUNDS_KEGG_IDS:
            if self._FindCompoundIndex(self.substrates, v_id) is not None:
                return True
            if self._FindCompoundIndex(self.products, v_id) is not None:
                return True
        return False

    def GetVolatileReactants(self):
        volatiles = []
        for v_id in VOLATILE_COMPOUNDS_KEGG_IDS:
            ind = self._FindCompoundIndex(self.substrates, v_id)
            if ind is not None:
                volatiles.append(self.substrates[ind])
            
        for v_id in VOLATILE_COMPOUNDS_KEGG_IDS:
            ind = self._FindCompoundIndex(self.products, v_id)
            if ind is not None:
                volatiles.append(self.products[ind])
        
        return volatiles
            
    def IsReactantFormulaMissing(self):
        try:
            self._GetAtomDiff()
            return False
        except ReactantFormulaMissingError:
            return True
    
    def IsBalanced(self):
        """Checks if the collection is atom-wise balanced.
        
        Returns:
            True if the collection is atom-wise balanced.
        """
        try:
            return self._IsBalanced(self._GetAtomDiff())
        except ReactantFormulaMissingError:
            return True
    
    def IsElectronBalanced(self):
        """Checks if the collection is electron-wise balanced.
        
        Returns:
            True if the collection is electron-wise balanced.
        """
        return self._GetElectronDiff() == 0
    
    def IsHalfReaction(self):
        """Checks if the reaction is a half-reaction (excess electrons).
        
        Returns:
            True if the collection is electron-wise balanced.
        """
        return self._GetElectronDiff() != 0
    
    def StandardizeHalfReaction(self):
        """Checks if the reaction is a half-reaction (excess electrons).
        
        Returns:
            True if the collection is electron-wise balanced.
        """
        if self._GetElectronDiff() > 0:
            self.SwapSides()
            
    def E0_tag(self):
        """Returns the standard transformed reduction potential of this reaction."""
        delta_electrons = abs(self._GetElectronDiff())
        assert delta_electrons != 0
        return - self.dg0_tag / (constants.F*delta_electrons)
    
    def _ExtraWaters(self):
        atom_diff = self._GetAtomDiff()
            
        if not atom_diff:
            return None
                
        # Ignore hydrogen.
        atom_diff.pop('H', 0)
        
        # Omit oxygen for checking balancedness.
        oxy_count = atom_diff.pop('O', 0)
    
        # If it's not balanced without oxygen, can't balance with water.
        if not self._IsBalanced(atom_diff):
            return None
        
        # Requires this many waters to balance (1 O per).
        return oxy_count

    @staticmethod
    def _FindCompoundIndex(side, id):
        """Returns the index of the compound with the given id.
        
        Args:
            side: a list of CompoundWithCoeff objects.
        
        Returns:
            The index of the compound or None if not present.
        """
        for i, c in enumerate(side):
            if c.compound.kegg_id == id:
                return i
        return None

    @staticmethod
    def _FindWater(side):
        """Returns the index of water into the list."""
        return Reaction._FindCompoundIndex(side, 'C00001')

    @staticmethod
    def _ReplaceCompound(side, from_id, to_id):
        """
            Replace a compound in the reaction with another compound according
            to their IDs.
            The stoichiometric coefficient and concentration are copied from
            the old compound to the new one.
        """
        index = Reaction._FindCompoundIndex(side, from_id)
        if index is None:
            return
        from_c_w_c = side[index]
        coeff = from_c_w_c.coeff
        conc = from_c_w_c.concentration
        to_c_w_c = CompoundWithCoeff.FromId(coeff, to_id, concentration=conc)
        side[index] = to_c_w_c

    @staticmethod
    def _AddCompound(side, id, how_many):
        """Adds "how_many" of the compound with the given id.
        
        Args:
            side: a list of CompoundWithCoeff objects.
            id: the KEGG id.
            how_many: how many waters to add.
        """
        i = Reaction._FindCompoundIndex(side, id)        
        if i:
            side[i].coeff += how_many
        else:
            c_w_coeff = CompoundWithCoeff.FromId(how_many, id)
            side.append(c_w_coeff)    

    @staticmethod
    def AddWater(side, how_many):
        """Adds "how_many" waters to a reaction side.
        
        Args:
            side: a list of CompoundWithCoeff objects.
            how_many: how many waters to add.
        """
        Reaction._AddCompound(side, 'C00001', how_many)
    
    @staticmethod
    def SubtractWater(side, how_many):
        """Removes at most "how_many" waters from a reaction side.
        
        Args:
            side: a list of CompoundWithCoeff objects.
            how_many: how many waters to subtract.
        
        Returns:
            How many waters are left after we subtracted as many as we could.
        """
        i = Reaction._FindWater(side)

        # Didn't find water in this side at all.        
        if i is None:
            return how_many
        
        net_water = side[i].coeff - how_many
        if net_water > 0:
            # We didn't consume all the waters on this side.
            side[i].coeff = net_water
            return 0
        
        # We consumed all the waters on this side.
        side.pop(i)
        return net_water
            
    def CanBalanceWithWater(self):
        """Returns True if balanced with or without water."""
        try:
            extra_waters = self._ExtraWaters()
        except ReactantFormulaMissingError as e:
            return True
        
        if extra_waters == None:
            return False
        return True

    @staticmethod
    def _FilterZeroes(side):
        """Removes compounds with coefficients equal to zero.
        
        Args:
            side: the reaction side to filter.
        """
        return filter(lambda x: x.coeff != 0, side)

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
        
        abs_waters = abs(extra_waters)
        if extra_waters > 0:
            waters_left = self.SubtractWater(self.substrates, abs_waters)
            self.AddWater(self.products, waters_left)
        else:
            waters_left = self.SubtractWater(self.products, abs_waters)
            self.AddWater(self.substrates, waters_left)
        
        self.substrates = self._FilterZeroes(self.substrates)
        self.products  = self._FilterZeroes(self.products)
        
        return True
    
    def TryReplaceCO2(self):
        """Attempt to replace CO2(aq) with CO2(total).
        
        Returns:
            True on success.
        """
        co2_id = 'C00011'
        bic_id = 'C00288'
        self._ReplaceCompound(self.substrates, co2_id, bic_id)
        self._ReplaceCompound(self.products, co2_id, bic_id)
        
        return self.TryBalanceWithWater()

    def BalanceElectrons(self,
                         acceptor_id='C00003',           # NAD+
                         reduced_acceptor_id='C00004'):  # NADH
        """Try to balance the reaction electons."""        
        net_electrons = self._GetElectronDiff()
        if net_electrons == 0:
            return
        
        acceptor = models.Compound.objects.get(kegg_id=acceptor_id)
        reduced_acceptor = models.Compound.objects.get(kegg_id=reduced_acceptor_id)
        
        if net_electrons < 0:
            # More product electrons. Need a donor on the left.
            num = (-net_electrons) / 2
            self._AddCompound(self.substrates, reduced_acceptor_id, num)
            self._AddCompound(self.products, acceptor_id, num)
        else:
            # More substrate-side electrons. Need an acceptor on the left.
            num = net_electrons / 2
            self._AddCompound(self.substrates, acceptor_id, num)
            self._AddCompound(self.products, reduced_acceptor_id, num)
        
    @staticmethod
    def _FilterHydrogen(compounds_with_coeffs):
        """Removes Hydrogens from the list of compounds."""
        return filter(lambda c: c.compound.kegg_id != 'C00080', compounds_with_coeffs)

    @staticmethod
    def GetTotalFormationEnergy(collection,
                                pH=constants.DEFAULT_PH,
                                pMg=constants.DEFAULT_PMG,
                                ionic_strength=constants.DEFAULT_IONIC_STRENGTH):
        """Compute an estimate for a collection of compounds + coefficients.
        
        You can compute the DeltaG of a reaction using negative coefficients for
        products.
        
        Args:
            collection: an iterable of CompoundWithCoeff objects.
            pH: the pH.
            pMg: the pMg.
            ionic_strength: the ionic strength.
        """ 
        sum = 0
        for compound_w_coeff in Reaction._FilterHydrogen(collection):
            c = compound_w_coeff.compound
            coeff = compound_w_coeff.coeff
            
            est = c.DeltaG(pH=pH, pMg=pMg, ionic_strength=ionic_strength)
            if est == None:
                logging.info('No estimate for compound %s', c.kegg_id)
                return None
            
            sum += coeff * est
    
        return sum
    
    def _GetConcentrationCorrection(self):
        """Get the concentration term in DeltaG' for these concentrations.
        
        Args:
            concentration_profile: a ConcentrationProfile object.
        
        Returns:
            The correction or None on error.
        """        
        # Ignore hydrogen for computing concentration corrections ala Alberty.
        rs = self._FilterHydrogen(self.substrates)
        ps = self._FilterHydrogen(self.products)
        
        # Shorthand for coeff * log(concentration)
        mult_log = lambda c: c.coeff * numpy.log(c.concentration)

        # Compute product and substrate terms.
        substrate_terms = [mult_log(c) for c in rs]
        product_terms = [mult_log(c) for c in ps]
        substrate_term = sum(substrate_terms)
        product_term = sum(product_terms)
        
        _r = constants.R
        _t = constants.DEFAULT_TEMP
        return _r * _t * (product_term - substrate_term)

    def DeltaG0(self):
        """Compute the DeltaG0 for a reaction.
        
        Returns:
            The DeltaG0 for this reaction, or None if data was missing.
        """
        substrates_sum = self.GetTotalFormationEnergy(
            self.substrates, pH=0, pMg=0, ionic_strength=0)
        products_sum = self.GetTotalFormationEnergy(
            self.products, pH=0, pMg=0, ionic_strength=0)
        if products_sum is None:
            logging.warning('Failed to get products formation energy.')
            return None
        if substrates_sum is None:
            print substrates_sum
            logging.warning('Failed to get substrates formation energy.')
            return None
        
        dg0 = products_sum - substrates_sum
        return dg0  

    def DeltaG0Tag(self, pH=None, pMg=None, ionic_strength=None):
        """Compute the DeltaG0' for a reaction.
        
        Returns:
            The DeltaG0' for this reaction, or None if data was missing.
        """
        ph = pH or self.ph
        pmg = pMg or self.pmg
        i_s = ionic_strength or self.i_s
        substrates_sum = self.GetTotalFormationEnergy(
            self.substrates, pH=ph, pMg=pmg, ionic_strength=i_s)
        products_sum = self.GetTotalFormationEnergy(
            self.products, pH=ph, pMg=pmg, ionic_strength=i_s)
        if products_sum is None:
            logging.warning('Failed to get products formation energy.')
            return None
        if substrates_sum is None:
            logging.warning('Failed to get substrates formation energy.')
            return None
        
        dg0_tag = products_sum - substrates_sum
        return dg0_tag

    def DeltaGTag(self, pH=None, pMg=None, ionic_strength=None):
        """Compute the DeltaG' for a reaction.
        
        Returns:
            The DeltaG' for this reaction, or None if data was missing.
        """
        dg0_tag = self.DeltaG0Tag(pH=pH, pMg=pMg, ionic_strength=ionic_strength)
        correction = self._GetConcentrationCorrection()
        return dg0_tag + correction
    
    def KeqTag(self):
        """Returns K'eq for this reaction."""
        dg0_tag = self.DeltaG0Tag()
        if dg0_tag is None:
            return None
        
        rt = constants.R * constants.DEFAULT_TEMP
        keq = numpy.exp(-dg0_tag / rt)
        return keq
    
    def NoDeltaGExplanation(self):
        """Get an explanation for why there's no delta G value.
        
        Return:
            The explanation or None.
        """
        for compound in self.substrates + self.products:
            if compound.compound.no_dg_explanation:
                name = compound.compound.common_names.all()[0].name
                return '%s %s' % (name,
                                  compound.compound.no_dg_explanation.lower())
        return None

    def AllCompoundsWithTransformedEnergies(self):
        for c_w_coeff in self.filtered_reactants:
            dgt = c_w_coeff.compound.DeltaG(pH=self.ph,
                                            pMg=self.pmg,
                                            ionic_strength=self.i_s)
            c_w_coeff.transformed_energy = dgt
            yield c_w_coeff

    def ExtraAtoms(self):
        try:
            diff = self._GetAtomDiff()
        except ReactantFormulaMissingError:
            return None
        diff.pop('H', 0)
        extras = filter(lambda t: t[1] > 0, diff.iteritems())
        if not extras:
            return None
        
        extras.sort(key=lambda t: t[1], reverse=True)
        return extras

    def MissingAtoms(self):
        try:
            diff = self._GetAtomDiff()
        except ReactantFormulaMissingError:
            return None
        diff.pop('H', 0)
        short = filter(lambda t: t[1] < 0, diff.iteritems())
        if not short:
            return None
        
        short = [(atom, -count) for atom, count in short]
        short.sort(key=lambda t: t[1], reverse=True)        
        return short
    
    def ExtraElectrons(self):
        diff = self._GetElectronDiff()
        if diff > 0:
            return diff
        return None
    
    def MissingElectrons(self):
        diff = self._GetElectronDiff()
        if diff < 0:
            return -diff
        return None
    
    def _CheckConservationLaw(self, sparse_reaction, claw):
        inner_prod = 0.0
        for kegg_id, coeff in claw.GetSparseRepresentation().iteritems():
            if kegg_id in sparse_reaction:
                inner_prod += coeff * sparse_reaction[kegg_id]
        return abs(inner_prod) < 1e-10
    
    def CheckConservationLaws(self):
        sparse_reaction = {}
        for compound in self.substrates:
            sparse_reaction[compound.compound.kegg_id] = -compound.coeff
        for compound in self.products:
            sparse_reaction[compound.compound.kegg_id] = compound.coeff

        all_claws = models.ConservationLaw.objects.select_related().all()
        for claw in all_claws:
            if not self._CheckConservationLaw(sparse_reaction, claw):
                return False
        return True
    
    contains_co2 = property(ContainsCO2)
    contains_volatile = property(ContainsVolatile)
    volatile_reactants = property(GetVolatileReactants)
    is_conserving = property(CheckConservationLaws)
    is_reactant_formula_missing = property(IsReactantFormulaMissing)
    is_balanced = property(IsBalanced)
    is_electron_balanced = property(IsElectronBalanced)
    is_half_reaction = property(IsHalfReaction)
    balanced_with_water = property(CanBalanceWithWater)
    extra_atoms = property(ExtraAtoms)
    missing_atoms = property(MissingAtoms)
    extra_electrons = property(ExtraElectrons)
    missing_electrons = property(MissingElectrons)
    filtered_reactants = property(lambda s: s._FilterHydrogen(s.substrates + s.products))
    filtered_substrates = property(lambda s: s._FilterHydrogen(s.substrates))
    filtered_products = property(lambda s: s._FilterHydrogen(s.products))
    all_compounds = property(AllCompoundsWithTransformedEnergies)
    dg0 = property(DeltaG0)
    dg0_tag = property(DeltaG0Tag)
    dg_tag = property(DeltaGTag)
    k_eq_tag = property(KeqTag)
    e0_tag = property(E0_tag)
    no_dg_explanation = property(NoDeltaGExplanation)
    standard_concentrations = property(StandardConcentrations)
    ph_graph_link = property(GetPhGraphLink)
    pmg_graph_link = property(GetPMgGraphLink)
    is_graph_link = property(GetIonicStrengthGraphLink)
    