#!/usr/bin/python

import json
import logging
import types
import re
import hashlib

import kegg_utils

from pygibbs import kegg_errors, kegg_parser


class Reaction(object):
    """A reaction from KEGG."""
    free_rid = -1 # class static variable
    
    def __init__(self, names, sparse_reaction,
                 rid=None, direction='<=>', weight=1):
        """Initialize the reaction."""
        self.SetNames(names)
        self.sparse = sparse_reaction
        if rid is None:
            self.rid = Reaction.free_rid
            Reaction.free_rid -= 1
        else:
            self.rid = rid
        self.weight = weight
        self.direction = direction
        self.definition = None
        self.equation = None
        self.ec_list = ['-.-.-.-']
    
    def __imul__(self, const):
        mult = float(const)
        for cid in self.sparse.iterkeys():
            self.sparse[cid] *= mult
        return self

    def __mul__(self, const):
        r = self.clone()
        r *= const
        return r
    
    def __rmul__(self, const):
        r = self.clone()
        r *= const
        return r
    
    def __iadd__(self, other):
        for cid in other.sparse.keys():
            if cid in self.sparse:
                self.sparse[cid] += other.sparse[cid]
            else:
                self.sparse[cid] = other.sparse[cid]
        
        zero_cids = [cid for cid in self.sparse.keys() if self.sparse[cid] == 0]
        for cid in zero_cids:
            del self.sparse[cid]
        return self
    
    def __add__(self, r1, r2):
        r = r1.clone()
        r += r2
        return r

    def SetNames(self, names):
        if type(names) == types.ListType:
            self.names = names
            self.name = names[0]
        elif type(names) in types.StringTypes:
            self.names = [names]
            self.name = names
        else:
            raise ValueError("Reaction name must be string or list of strings")

    def clone(self):
        reaction = Reaction(self.names, dict(self.sparse))
        reaction.rid = self.rid
        reaction.weight = self.weight
        reaction.direction = self.direction
        reaction.definition = self.definition
        reaction.equation = self.equation
        reaction.ec_list = self.ec_list
        return reaction

    def reverse(self):
        """ Returns the reverse Reaction with, i.e. substrates become products
            and vice versa.
        """
        reaction = self.clone()
        reaction.sparse = dict([(cid, -coeff) for 
                                (cid, coeff) in self.iteritems()])
        if self.direction == '<=>':
            reaction.direction = '<=>'
        elif self.direction == '=>':
            reaction.direction = '<='
        elif self.direction == '<=':
            reaction.direction = '=>'
        else:
            raise Exception('invalid direction: ' + self.direction)
        return reaction
    
    def Normalize(self):
        """
            Set the coefficient of the reactant with the lowest ID to 1,
            by multiplying all coefficients by its reciprocal
        """
        factor = 1.0 / self.sparse[min(self.sparse.keys())]
        self.sparse = dict((k,v*factor) for (k,v) in self.iteritems())
        return factor
    
    @staticmethod
    def FromEntryDict(key, field_map):
        if not key.startswith('R'):
            return None
            
        equation_value = field_map.GetStringField("EQUATION", default_value="<=>")
        try:
            reaction = Reaction.FromFormula(equation_value)
            reaction.SetNames(key)
            reaction.rid = int(key[1:])

            if "ENZYME" in field_map:
                reaction.ec_list = field_map.GetStringListField('ENZYME')
            if "NAME" in field_map:
                reaction.SetNames(kegg_parser.NormalizeNames(field_map["NAME"]))
            if "DEFINITION" in field_map:
                reaction.definition = field_map["DEFINITION"]
            if "EQUATION" in field_map:
                reaction.equation = field_map["EQUATION"]
        except (kegg_errors.KeggParseException,
                kegg_errors.KeggNonCompoundException,
                ValueError):
            logging.debug("Cannot parse reaction formula: " + equation_value)
            return None
        
        return reaction
    
    def ToDBRow(self, rid):
        """Build a database row from a reaction."""
        return [rid,
                json.dumps(self.names),
                self.definition,
                json.dumps(self.ec_list),
                self.equation]
        
    @staticmethod
    def FromDBRow(row_dict):
        """Build a Reaction from a database row."""
        reaction = Reaction.FromFormula(row_dict['equation'])
        names = json.loads(row_dict['all_names'])
        reaction.SetNames(names)
        reaction.rid = row_dict['rid']
        reaction.equation = row_dict['equation']
        reaction.definition = row_dict['definition']
        reaction.ec_list = json.loads(row_dict['ec_list'])
        return reaction
    
    @staticmethod
    def parse_reaction_formula_side(s):
        """ 
            Parses the side formula, e.g. '2 C00001 + C00002 + 3 C00003'
            Ignores stoichiometry.
            
            Returns:
                The set of CIDs.
        """
        if s.strip() == "null":
            return {}
        
        compound_bag = {}
        for member in re.split('\s+\+\s+', s):
            tokens = member.split(None, 1)
            if len(tokens) == 1:
                amount = 1
                key = member
            else:
                try:
                    amount = float(tokens[0])
                except ValueError:
                    raise kegg_errors.KeggParseException(
                        "Non-specific reaction: %s" % s)
                key = tokens[1]
                
            if key[0] != 'C':
                raise kegg_errors.KeggNonCompoundException(
                    "Compound ID doesn't start with C: %s" % key)
            try:
                cid = int(key[1:])
                compound_bag[cid] = compound_bag.get(cid, 0) + amount
            except ValueError:
                raise kegg_errors.KeggParseException(
                    "Non-specific reaction: %s" % s)
        
        return compound_bag

    @staticmethod
    def FromFormula(formula, name='reaction', rid=None):
        """ 
            Parses a two-sided formula such as: 2 C00001 => C00002 + C00003 
            
            Return:
                The set of substrates, products and the direction of the reaction
        """
        tokens = re.findall("([^=^<]+) (<*=>*) ([^=^>]+)", formula)
        if len(tokens) != 1:
            raise kegg_errors.KeggParseException(
                "Cannot parse this formula: %s" % formula)
        
        left, direction, right = tokens[0] # the direction: <=, => or <=>
        
        sparse_reaction = {}
        for cid, count in Reaction.parse_reaction_formula_side(left).iteritems():
            sparse_reaction[cid] = sparse_reaction.get(cid, 0) - count 
    
        for cid, count in Reaction.parse_reaction_formula_side(right).iteritems():
            sparse_reaction[cid] = sparse_reaction.get(cid, 0) + count 
    
        reaction = Reaction(name, sparse_reaction, rid=rid, direction=direction)
        return reaction

    def replace_compound(self, replace_cid, with_cid):
        """Replace one CID with another in this reaction.
        
        Args:
            replace_cid: the CID to replace.
            with_cid: the one to replace with.
        """
        if replace_cid not in self.sparse:
            return
        
        if with_cid in self.sparse:
            raise ValueError('Reaction %s already contains CID %s' % (self.rid,
                                                                      with_cid))
        
        count = self.sparse.pop(replace_cid)
        self.sparse[with_cid] = count

    def get_cids(self):
        """Returns the KEGG IDs of the products and reactants."""
        return set(self.sparse.keys())

    @staticmethod
    def BalanceSparseReaction(sparse, balance_water=False,
                              balance_hydrogens=False, exception_if_unknown=False):
        from pygibbs.kegg import Kegg
        kegg = Kegg.getInstance()
        kegg_utils.balance_reaction(kegg, sparse, balance_water,
                                    balance_hydrogens, exception_if_unknown)

    def Balance(self, balance_water=False, balance_hydrogens=False,
                exception_if_unknown=False):
        """
            Balances a reaction
            
            Arguments:
                If balance_water=True and there is an imbalance of oxygen atoms, Balance
                changes the reaction by adding H2O until it is balanced.
                If balance_hydrogens=True then H+ are used to balance the amount of hydrogen atoms.
                
                If exception_if_unknown=True then an exception will be raised also if there
                    is not enough information to know if the reaction is balanced or not.
            
            If the reaction cannot be balanced, raises KeggReactionNotBalancedException
        """
        Reaction.BalanceSparseReaction(self.sparse, balance_water, 
                                       balance_hydrogens, exception_if_unknown)
        
    def PredictReactionEnergy(self, thermodynamics, 
                              pH=None, pMg=None, I=None ,T=None):
        return float(thermodynamics.GetTransfromedKeggReactionEnergies([self], pH=pH, pMg=pMg, I=I, T=T))
    
    def HashableReactionString(self):
        """
            Return a hashable string for a biochemical reaction.
        
            The string fully identifies the biochemical reaction including its direction.
            If it is equal to another reaction's string, then they have identical
            stoichiometry and direction.
        """
        sort_key = lambda r: r[0]
        make_str = lambda r: '%d %.2f' % r
        is_not_hydrogen = lambda r: r[0] != 80
        
        reactants_strs = map(make_str,
                             sorted(filter(is_not_hydrogen, self.iteritems()),
                                    key=sort_key))
        return ' + '.join(reactants_strs)
    
    @staticmethod
    def HashReaction(reaction):
        md5 = hashlib.md5()
        md5.update(reaction.HashableReactionString())
        return md5.hexdigest()
    
    def __hash__(self):
        return hash(Reaction.HashReaction(self))
    
    def __eq__(self, other):
        return self.HashableReactionString() == other.HashableReactionString()
    
    @staticmethod
    def write_compound_and_coeff(cid, coeff, show_cids=True):
        if show_cids:
            comp = "C%05d" % cid
        else:
            from pygibbs.kegg import Kegg
            kegg = Kegg.getInstance()
            comp = kegg.cid2name(cid)
        if coeff == 1:
            return comp
        else:
            return "%g %s" % (coeff, comp)

    @staticmethod
    def write_full_reaction(sparse, direction='=>', show_cids=True):
        """String representation."""
        left = []
        right = []
        for cid, coeff in sorted(sparse.iteritems()):
            if coeff < 0:
                left.append(Reaction.write_compound_and_coeff(cid, -coeff, show_cids=show_cids))
            elif coeff > 0:
                right.append(Reaction.write_compound_and_coeff(cid, coeff, show_cids=show_cids))
        return "%s %s %s" % (' + '.join(left), direction, ' + '.join(right))

    def FullReactionString(self, show_cids=True):
        return self.write_full_reaction(self.sparse, self.direction, show_cids=show_cids)

    def __str__(self):
        return self.name + ': ' + self.FullReactionString()
    
    def iteritems(self):
        return self.sparse.iteritems()
    
    def to_hypertext(self, show_cids=True):
        from pygibbs.kegg import Kegg
        my_kegg = Kegg.getInstance()
        return my_kegg.sparse_to_hypertext(self.sparse, show_cids=show_cids)
    def to_hypertext_with_names(self):
        return self.to_hypertext(show_cids=False)
    hypertext = property(to_hypertext)
    hypertext_with_names = property(to_hypertext_with_names)
    
    def is_not_futile(self):
        return max([abs(x) for x in self.sparse.values()]) > 0.01
    
    def get_link(self):
        return ('http://www.genome.jp/dbget-bin/www_bget?rn:R%05d'
                % self.rid)
        
    def ToJSONDict(self):
        """Format the reaction as a JSON dictionary."""
        reaction = [(coeff, 'C%05d' % cid) for cid, coeff
                    in sorted(self.iteritems())]
        d = {'names': self.names, 'reaction': reaction}
        d['RID'] = 'R%05d' % self.rid
        if self.ec_list:
            d['ECS'] = self.ec_list
        if self.equation:
            d['equation'] = self.equation
        if self.weight:
            d['weight'] = self.weight
        if self.direction:
            d['direction'] = self.direction
        if self.definition:
            d['definition'] = self.definition
            
        return d

    def ToJSONString(self):
        return json.dumps(self.ToJSONDict())
    
    @staticmethod
    def FromJSONString(s):
        d = json.loads(s)
        sparse = {}
        for coeff, cid in d['reaction']:
            sparse[int(cid[1:])] = float(coeff)
        rid = int(d['RID'][1:])
        reaction = Reaction(d['names'], sparse, rid=rid)
        if 'ECS' in d:
            reaction.ec_list = d['ECS']
        if 'equation' in d:
            reaction.equation = d['equation']
        if 'weight' in d:
            reaction.weight = d['weight']
        if 'direction' in d:
            reaction.direction = d['direction']
        if 'definition' in d:
            reaction.definition = d['definition']
        return reaction

def GetAllReactionsFromDB(db):
    """Fetch all the compounds from the database."""
    reaction_list = []
    for row_dict in db.DictReader('kegg_reaction'):
        reaction_list.append(Reaction.FromDBRow(row_dict))
    return reaction_list
