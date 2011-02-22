#!/usr/bin/python


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

    def __str__(self):
        """String representation."""
        left = []
        right = []
        for (cid, coeff) in sorted(self.sparse.iteritems()):
            if (coeff < 0):
                left.append(self.write_compound_and_coeff(cid, -coeff))
            elif (coeff > 0):
                right.append(self.write_compound_and_coeff(cid, coeff))
        return "%s -> %s" % (' + '.join(left), ' + '.join(right))
    
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