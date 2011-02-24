#!/usr/bin/python
from pygibbs import kegg_utils


class Reaction(object):
    """A reaction from KEGG."""
    free_rid = -1 # class static variable
    
    def __init__(self, names, sparse_reaction,
                 rid=None, direction='<=>', weight=1):
        """Initialize the reaction."""
        self.names = names
        self.name = names[0]
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

    def get_cids(self):
        """Returns the KEGG IDs of the products and reactants."""
        return set(self.sparse.keys())
    
    def unique_string(self):
        """Returns a unique string for the reaction."""
        formatted = ["%d C%05d" % (coeff, cid) for (cid, coeff)
                     in sorted(self.sparse.iteritems())]
        return " + ".join(formatted)

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
        for (cid, coeff) in sorted(sparse.iteritems()):
            if (coeff < 0):
                left.append(Reaction.write_compound_and_coeff(cid, -coeff))
            elif (coeff > 0):
                right.append(Reaction.write_compound_and_coeff(cid, coeff))
        return "%s -> %s" % (' + '.join(left), ' + '.join(right))
        

    def __str__(self):
        return Reaction.write_full_reaction(self.sparse)
    
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
