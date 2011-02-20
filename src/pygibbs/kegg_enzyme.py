#!/usr/bin/python

from pygibbs import kegg_parser

class Enzyme(object):
    """A class representing an Enzyme from Kegg."""
    
    def __init__(self, ec_class, title=None, names=None,
                 reactions=None, substrate=None, product=None,
                 cofactor=None, organism=None):
        """Initialize the enzyme class.
        
        Args:
            ec_class: the EC class of the enzyme. Used as the primary identifier.
            title: the title field from KEGG.
            names: common names for the enzyme.
            reactions: a list of (integer) KEGG reaction IDs.
            substrate: the substrate of the reaction, if defined.
            product: the product of the enzyme, if defined.
            cofactor: cofactors used by the enzyme.
            organism: the organism this entry is found in, if defined.
            
        Attributes:
            ec
            title
            names
            reactions
            substrate
            product
            cofactor
            organism
        """
        self.ec = ec_class
        self.title = title
        self.names = names or []
        self.reactions = reactions or []
        self.substrate = substrate
        self.product = product
        self.cofactor = cofactor
        self.organism = organism

    @staticmethod
    def ProcessEC(ec_str):
        if ec_str.startswith('EC '):
            return ec_str[3:]
        return ec_str

    @staticmethod
    def FromEntryDict(ec_class, entry_dict):
        """Initialize from an kegg_parser.EntryDictWrapper instance.
        
        Args:
            ec_class: the EC class of the enzyme. Used as the primary identifier.
            entry_dict. A kegg_parser.EntryDictWrapper instance.
            
        Returns:
            A new Enzyme object with appropriately initialized parameters.
        """
        enz = Enzyme(Enzyme.ProcessEC(ec_class))
        if 'NAME' in entry_dict:
            enz.names = kegg_parser.NormalizeNames(entry_dict.get("NAME"))
        if 'ALL_REAC' in entry_dict:
            enz.reactions = kegg_parser.NormalizeReactions(
                entry_dict.get('ALL_REAC'))
        enz.substrate = entry_dict.get('SUBSTRATE', None)
        enz.product = entry_dict.get('PRODUCT', None)
        enz.title = entry_dict.get('TITLE', None)
        enz.cofactor = entry_dict.get('COFACTOR', None)
        enz.organism = entry_dict.get('ORGANISM', None)
        return enz

    @staticmethod
    def FromDBRow(row):
        """Initialize from a dictionary representing a DB row.
        
        The row was (hopefully) created by calling ToDBRow on an Enzyme instance.
        
        Args:
            row: The row read from the DB (as a dictionary).
            
        Returns:
            A new Enzyme object with appropriately initialized parameters.
        """
        if not 'ec' in row:
            return None
        
        enz = Enzyme(Enzyme.ProcessEC(row['ec']))
        if 'all_names' in row:
            names = row['all_names']
            if names:
                enz.names = names.split(', ')
        if 'rid_list' in row:
            rid_list = row['rid_list']
            if rid_list:
                enz.reactions = [int(r) for r in row['rid_list'].split(', ')]
        enz.substrate = row.get('substrate', None)
        enz.product = row.get('product', None)
        enz.title = row.get('title', None)
        enz.cofactor = row.get('cofactor', None)
        enz.organism = row.get('organism', None)
        return enz

    def ToDBRow(self):
        """Create a DB row from this enzyme."""
        return [self.ec,
                ', '.join(self.names),
                self.title,
                ', '.join([str(r) for r in self.reactions]),
                self.substrate,
                self.product,
                self.cofactor,
                self.organism]

    def __str__(self):
        """String representation of the enzyme."""
        l = []
        if self.ec:
            l.append('EC: %s\n' % self.ec)
        if self.title:
            l.append('Title: %s\n' % self.title)
        if self.names:
            l.append('Names: %s\n' % ', '.join(self.names))
        if self.reactions:
            l.append('Reactions: %s\n' % ', '.join(self.reactions))
        if self.substrate:
            l.append('Substrate: %s\n' % self.substrate)
        if self.product:
            l.append('Product: %s\n' % self.product)
        if self.cofactor:
            l.append('Co-factor: %s\n' % self.cofactor)
        if self.organism:
            l.append('Organism: %s\n' % self.organism)
        
        return ''.join(l)
    
    def GetLink(self):
        """Returns a link to the KEGG page for this enzyme."""
        return 'http://kegg.jp/dbget-bin/www_bget?ec:%s' % self.ec
    
    def HasReactions(self):
        """Returns True if this enzyme has reactions defined."""
        if self.reactions:
            return True
        return False