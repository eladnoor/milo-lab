#!/usr/bin/python

import re
import json

from pygibbs import kegg_parser


class Gene(object):
    """A class representing a Gene from Kegg."""
    
    def __init__(self, organism, gene_id,
                 aa_size=None, nt_size=None):
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
        self.organism = organism
        self.gene_id = gene_id
        self.aa_size = aa_size
        self.nt_size = nt_size
        
    @staticmethod
    def FromEntryDict(entry_dict):
        """Initialize from an kegg_parser.EntryDictWrapper instance.
        
        Args:
            ec_class: the EC class of the enzyme. Used as the primary identifier.
            entry_dict. A kegg_parser.EntryDictWrapper instance.
            
        Returns:
            A new Enzyme object with appropriately initialized parameters.
        """
        organism = entry_dict.get('ORGANISM')
        gene_id = entry_dict.get('ENTRY').split()[0]
        aa_size = entry_dict.get('AASEQ').split()[0]
        nt_size = entry_dict.get('NTSEQ').split()[0]
        aa_size = int(aa_size)
        nt_size = int(nt_size)
        
        return Gene(organism, gene_id, aa_size=aa_size, nt_size=nt_size)
        
    @staticmethod
    def FromDBRow(row_dict):
        """Initialize from a dictionary representing a DB row_dict.
        
        The row_dict was (hopefully) created by calling ToDBRow on an Enzyme instance.
        
        Args:
            row_dict: The row_dict read from the DB (as a dictionary).
            
        Returns:
            A new Enzyme object with appropriately initialized parameters.
        """
        # TODO(flamholz): implement
        raise NotImplementedError

    def ToDBRow(self):
        """Create a DB row from this enzyme."""
        # TODO(flamholz): implement
        raise NotImplementedError

    def __str__(self):
        """String representation of the enzyme."""
        l = []
        if self.organism:
            l.append('Organism: %s\n' % self.organism)
        if self.gene_id:
            l.append('Gene ID: %s\n' % self.organism)
        if self.aa_size:
            l.append('AA Size: %d\n' % self.aa_size)
        if self.nt_size:
            l.append('NT Size: %d\n' % self.nt_size)

        return ''.join(l)
    