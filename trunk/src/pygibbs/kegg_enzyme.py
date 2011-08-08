#!/usr/bin/python

import re
import json

from pygibbs import kegg_parser


class Enzyme(object):
    """A class representing an Enzyme from Kegg."""
    
    def __init__(self, ec_class, title=None, names=None,
                 reactions=None, substrates=None, products=None,
                 cofactors=None, organisms=None, orthology=None,
                 genes=None):
        """Initialize the enzyme class.
        
        Args:
            ec_class: the EC class of the enzyme. Used as the primary identifier.
            title: the title field from KEGG.
            names: common names for the enzyme.
            reactions: a list of (integer) KEGG reaction IDs.
            substrates: the substrate of the reaction, if defined.
            products: the products of the enzyme, if defined.
            cofactors: cofactors used by the enzyme.
            organisms: the organisms this enzyme is found in, if defined.
            orthology: the mapping dictionary from orthology IDs to names.
            genes: the mapping from organisms to a list of gene names.
            
        Attributes:
            ec
            title
            names
            reactions
            substrates
            products
            cofactors
            organisms
            orthology
            genes
        """
        self.ec = ec_class
        self.title = title
        self.names = names or []
        self.reactions = reactions or []
        self.substrates = substrates 
        self.products = products
        self.cofactors = cofactors
        self.organisms = organisms or []
        self.orthology = orthology or {}
        self.genes = genes or {}

    @staticmethod
    def ProcessEC(ec_str):
        """Format EC strings properly.
        
        Args:
            ec_str: the raw EC string from the KEGG file. May include
              the "EC " prefix.
        """
        stripped = ec_str.strip()
        if stripped.startswith('EC '):
            return stripped[3:]
        return stripped

    @staticmethod
    def GetCompoundIds(cpd_str):
        """Returns a list of KEGG ids parsed from the string."""
        if not cpd_str:
            return None
        
        pattern = re.compile(r'\[CPD:(C\d{5})\]')
        return pattern.findall(cpd_str)
        
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
        if 'ORGANISM' in entry_dict:
            enz.organisms = kegg_parser.NormalizeOrganisms(
                entry_dict.get('ORGANISM'))
        if 'ORTHOLOGY' in entry_dict:
            enz.orthology = kegg_parser.ParseOrthologyMapping(
                entry_dict.get('ORTHOLOGY'))
        if 'GENES' in entry_dict:
            enz.genes = kegg_parser.ParseOrganismToGeneMapping(
                entry_dict.get('GENES'))
        
        enz.substrates = Enzyme.GetCompoundIds(entry_dict.get('SUBSTRATE', None))
        enz.products = Enzyme.GetCompoundIds(entry_dict.get('PRODUCT', None))
        enz.cofactors = Enzyme.GetCompoundIds(entry_dict.get('COFACTOR', None))
        enz.title = entry_dict.get('TITLE', None)
        
        if enz.title:
            enz.title = enz.title.replace('\t', ' ')

        return enz

    @staticmethod
    def FromDBRow(row_dict):
        """Initialize from a dictionary representing a DB row_dict.
        
        The row_dict was (hopefully) created by calling ToDBRow on an Enzyme instance.
        
        Args:
            row_dict: The row_dict read from the DB (as a dictionary).
            
        Returns:
            A new Enzyme object with appropriately initialized parameters.
        """
        if not 'ec' in row_dict:
            return None
        
        enz = Enzyme(Enzyme.ProcessEC(row_dict['ec']))
        if 'all_names' in row_dict:
            names = row_dict['all_names']
            if names:
                enz.names = json.loads(names)
                
        if 'rid_list' in row_dict:
            rid_list = row_dict['rid_list']
            if rid_list:
                enz.reactions = map(int, json.loads(rid_list))
                
        if 'organism_list' in row_dict:
            org_list = row_dict['organism_list']
            if org_list:
                enz.organisms = json.loads(org_list)
                        
        enz.substrates = row_dict.get('substrate', None)
        if enz.substrates:
            enz.substrates = json.loads(enz.substrates)
            
        enz.products = row_dict.get('product', None)
        if enz.products:
            enz.products = json.loads(enz.products)
            
        enz.cofactors = row_dict.get('cofactor', None)
        if enz.cofactors:
            enz.cofactors = json.loads(enz.cofactors)
            
        enz.title = row_dict.get('title', None)
        
        orthology = row_dict.get('orthology_map')
        genes = row_dict.get('genes_map')
        if orthology:
            enz.orthology = json.loads(orthology)
        if genes:
            enz.genes = json.loads(genes)
        
        return enz

    def ToDBRow(self):
        """Create a DB row from this enzyme."""
        row = [self.ec,
               json.dumps(self.names),
               self.title,
               json.dumps([str(r) for r in self.reactions])]
        
        if self.substrates:
            row.append(json.dumps(self.substrates))
        else:
            row.append(None)
        
        if self.products:
            row.append(json.dumps(self.products))
        else:
            row.append(None)
        
        if self.cofactors:
            row.append(json.dumps(self.cofactors))
        else:
            row.append(None)
            
        if self.organisms:
            row.append(json.dumps(self.organisms))
        else:
            row.append(None)
            
        if self.orthology:
            row.append(json.dumps(self.orthology))
        else:
            row.append(None)
            
        if self.genes:
            row.append(json.dumps(self.genes))
        else:
            row.append(None)
        
        return row

    @staticmethod
    def GetStringRID(int_rid):
        """Gets a string RID from an integer one."""
        return 'R%05d' % int_rid

    def ToJSONDict(self):
        """Format the enzyme as a JSON dictionary."""
        rids = map(self.GetStringRID, self.reactions)
        return {'EC': self.ec,
                'title': self.title,
                'names': self.names,
                'reaction_ids': rids,
                'substrates': self.substrates,
                'products': self.products,
                'cofactors': self.cofactors,
                'organisms': self.organisms,
                'orthology': self.orthology,
                'genes': self.genes}

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
        if self.substrates:
            l.append('Substrates: %s\n' % ', '.join(self.substrates))
        if self.products:
            l.append('Product: %s\n' % ', '.join(self.products))
        if self.cofactors:
            l.append('Co-factor: %s\n' % ', '.join(self.cofactors))
        if self.organisms:
            l.append('Organism: %s\n' % self.organisms)
        
        return ''.join(l)
    
    def GetLink(self):
        """Returns a link to the KEGG page for this enzyme."""
        return 'http://kegg.jp/dbget-bin/www_bget?ec:%s' % self.ec
    kegg_link = property(GetLink)
    
    def HasReactions(self):
        """Returns True if this enzyme has reactions defined."""
        if self.reactions:
            return True
        return False