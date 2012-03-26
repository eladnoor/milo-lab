#!/usr/bin/python

import sys
import numpy as np

from pygibbs import kegg
from pygibbs import kegg_parser
from pygibbs import kegg_enzyme
from pygibbs import kegg_gene
from SOAPpy import WSDL
from optparse import OptionParser


def MakeOpts():
    """Returns an OptionParser object with all the default options."""
    opt_parser = OptionParser()
    opt_parser.add_option("-e", "--enzyme_ec", 
                          dest="enzyme_ec",
                          help="The EC number of the enzyme")
    return opt_parser


def Main():
    serv = WSDL.Proxy(kegg.Kegg.WSDL_URL)
    
    options, _ = MakeOpts().parse_args(sys.argv)
    assert options.enzyme_ec
    
    key = 'EC %s' % options.enzyme_ec
    enz_data = serv.bget(key)
    entry2fields_map = kegg_parser.ParsedKeggFile.FromKeggAPI(enz_data)
    field_map = entry2fields_map.get(key)
    enz = kegg_enzyme.Enzyme.FromEntryDict(key, field_map)

    gene_keys = []
    for org, gene_ids in enz.genes.iteritems():
        gene_keys.extend(['%s:%s' % (org, gid) for gid in gene_ids])
    
    genes_query = ' '.join(gene_keys)
    genes_data = serv.bget(genes_query)
    entry2fields_map = kegg_parser.ParsedKeggFile.FromKeggAPI(genes_data)
    genes = []
    for key in sorted(entry2fields_map.keys()):
        field_map = entry2fields_map[key]
        genes.append(kegg_gene.Gene.FromEntryDict(field_map))

    aa_sizes = np.array([g.aa_size for g in genes])
    enz_name = '(%s)' % enz.ec
    if enz.names:
        enz_name = '%s %s' % (enz.names[0], enz_name)
    mean_size = np.mean(aa_sizes)
    geom_mean_size = np.exp(np.mean(np.log(aa_sizes)))
    print 'Enzyme:', enz_name
    print 'Geometric Mean AA Size', geom_mean_size
    print 'Average AA Size', mean_size
    print 'Std Dev', np.std(aa_sizes)
    

if __name__ == '__main__':
    Main()


