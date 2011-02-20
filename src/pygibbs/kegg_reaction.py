#!/usr/bin/python


class Reaction(object):
    """A reaction from KEGG."""
    free_rid = -1 # class static variable
    
    def __init__(self, name, sparse_reaction,
                 rid=None, direction='<=>', weight=1):
        """Initialize the reaction."""
        self.name = name
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
    
    def is_specific(self, cid2atom_bag):
        for cid in self.sparse.keys():
            atom_bag = cid2atom_bag.get(cid, None)
            if atom_bag == None or 'R' in atom_bag:
                # This reaction cannot be checked since there is an
                # unspecific compound
                return False
        return True
    
    def is_balanced(self, cid2atom_bag, balance_water=True):
        """
            Checks if the reaction is balanced: i.e. the sum of elements is conserved (not including hydrogen atoms).
            If oxygen is not balanced, the method adds water molecules in the right amount to fix it.
        """
        atom_diff = {}
        for cid, coeff in self.sparse.iteritems():
            atom_bag = cid2atom_bag.get(cid, None)
            for atomic_number, atom_count in atom_bag.iteritems():
                new_count =  atom_diff.get(atomic_number, 0)
                new_count += coeff * atom_count
                atom_diff[atomic_number] = new_count

        # ignore H and O inconsistencies
        if 'H' in atom_diff:
            del atom_diff['H']
        if balance_water:
            if 'O' in atom_diff and atom_diff['O'] != 0:
                self.sparse[1] = self.sparse.get(1, 0) - atom_diff['O']
                del atom_diff['O']
        
        return max([abs(x) for x in atom_diff.values()]) < 0.01
    
    def verify(self, cid2atom_bag):
        if not self.is_specific(cid2atom_bag):
            # unspecific reactions cannot be automatically checked for
            # chemical balance therefore we cannot use them here.
            return "unspecific"
        elif not self.is_balanced(cid2atom_bag):
            return "unbalanced"
        elif not self.is_not_futile():
            return "futile"
        else:
            return None
    
    def get_link(self):
        return ('http://www.genome.jp/dbget-bin/www_bget?rn:R%05d'
                % self.rid)