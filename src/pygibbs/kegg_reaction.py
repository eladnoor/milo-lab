#!/usr/bin/python
import kegg_utils
import hashlib
import types

class Reaction(object):
    """A reaction from KEGG."""
    free_rid = -1 # class static variable
    
    def __init__(self, names, sparse_reaction,
                 rid=None, direction='<=>', weight=1):
        """Initialize the reaction."""
        if type(names) == types.ListType:
            self.names = names
            self.name = names[0]
        elif type(names) == types.StringType:
            self.names = [names]
            self.name = names
        self.sparse = sparse_reaction
        if rid == None:
            self.rid = Reaction.free_rid
            Reaction.free_rid -= 1
        else:
            self.rid = rid
        self.weight = weight
        self.direction = direction
        self.definition = None
        self.equation = None
        self.ec_list = '-.-.-.-'

    def clone(self):
        reaction = Reaction(self.names, self.sparse)
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
        elif self.direction == '=>':
            reaction.direction = '=>'
        else:
            raise Exception('invalid direction: ' + self.direction)
        return reaction
        
    @staticmethod
    def FromDBRow(row_dict):
        """Build a Reaction from a database row."""
        (sparse, direction) = kegg_utils.parse_reaction_formula(row_dict['equation'])
        names = row_dict['all_names'].split(';')
        reaction = Reaction(names=names, sparse_reaction=sparse, 
                            rid=row_dict['rid'], direction=direction)
        reaction.equation = row_dict['equation']
        reaction.definition = row_dict['definition']
        reaction.ec_list = row_dict['ec_list']
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

    def Balance(self, balance_water=False):
        from pygibbs.kegg import Kegg
        kegg = Kegg.getInstance()
        self.sparse = kegg.BalanceReaction(self.sparse, balance_water)

    def PredictReactionEnergy(self, thermodynamics, 
                              pH=None, pMg=None, I=None ,T=None):
        return sum([coeff * thermodynamics.cid2dG0_tag(cid, pH, pMg, I, T)
                    for cid, coeff in self.sparse.iteritems()])
    
    def HashableReactionString(self):
        """
            Return a hashable string for a biochemical reaction.
        
            The string fully identifies the biochemical reaction including its direction.
            If it is equal to another reaction's string, then they have identical
            stoichiometry and direction.
        """
        sort_key = lambda r: r[0]
        make_str = lambda r: '%d %.2f' % r
        is_not_hydrogen = lambda r: r[0] != 'C00080'
        
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
    def write_compound_and_coeff(cid, coeff):
        if (coeff == 1):
            return "C%05d" % cid
        else:
            return "%g C%05d" % (coeff, cid)

    @staticmethod
    def write_full_reaction(sparse):
        """String representation."""
        left = []
        right = []
        for cid, coeff in sorted(sparse.iteritems()):
            if (coeff < 0):
                left.append(Reaction.write_compound_and_coeff(cid, -coeff))
            elif (coeff > 0):
                right.append(Reaction.write_compound_and_coeff(cid, coeff))
        return "%s -> %s" % (' + '.join(left), ' + '.join(right))

    def __str__(self):
        return self.name + ': ' + Reaction.write_full_reaction(self.sparse)
    
    def iteritems(self):
        return self.sparse.iteritems()
    
    def to_hypertext(self, show_cids=True):
        from pygibbs.kegg import Kegg
        kegg = Kegg.getInstance()
        return kegg.sparse_to_hypertext(self.sparse, show_cids=show_cids)
    
    def is_not_futile(self):
        return max([abs(x) for x in self.sparse.values()]) > 0.01
    
    def get_link(self):
        return ('http://www.genome.jp/dbget-bin/www_bget?rn:R%05d'
                % self.rid)
        
    def ToJSONDict(self):
        """Format the reaction as a JSON dictionary."""
        reaction = [(coeff, 'C%05d' % cid) for cid, coeff
                    in sorted(self.sparse.iteritems())]
        return {'RID': 'R%05d' % self.rid,
                'names': self.names,
                'ECS': self.ec_list,
                'reaction': reaction}

def GetAllReactionsFromDB(db):
    """Fetch all the compounds from the database."""
    reaction_list = []
    for row_dict in db.DictReader('kegg_reaction'):
        reaction_list.append(Reaction.FromDBRow(row_dict))
    return reaction_list
