#!/usr/bin/python

import sys
from collections import Counter
from optparse import OptionParser

from genomics import genome_db
from genomics import pathway

    
def MakeOpts():
    """Returns an OptionParser object with all the default options."""
    opt_parser = OptionParser()
    opt_parser.add_option("-p", "--pathways_filename",
                          dest="pathways_filename",
                          default='../data/genomics/glycolysis_pathways_unique.json',
                          help="Input pathways JSON file")
    opt_parser.add_option("-d", "--genome_db_filename",
                          dest="genome_db_filename",
                          default='../res/genomes.sqlite',
                          help="Genome database filename")
    return opt_parser


def Main():
    options, _ = MakeOpts().parse_args(sys.argv)
    assert options.pathways_filename and options.genome_db_filename
    
    pathways = pathway.LoadPathways(options.pathways_filename)
    db = genome_db.GenomeDB(options.genome_db_filename)
    
    for path in pathways:
        org_counts = {}
        for enz_set in path.enzyme_sets:
            orgs_w_enz = set()
            for ec in enz_set:
                q = db.db.Execute("SELECT organism FROM organism_enzymes WHERE EC='%s'" % ec)
                orgs_w_enz.update([i[0] for i in q])
                
            for org in orgs_w_enz:
                org_counts[org] = org_counts.get(org, 0) + 1
        
        orgs = [o for o,c in org_counts.iteritems()
                if c == len(path.enzyme_sets)]
        broad_oxygen_reqs = []
        energy_srcs = []
        metabolism = []
        for org in orgs:
            broad_oxygen_reqs.append(db.KEGG2BroadOxygenReq(org))
            energy_srcs.append(db.KEGG2EnergySource(org))
            metabolism.append(db.KEGG2Metabolism(org))
        req_counts = Counter(broad_oxygen_reqs)
        energy_counts = Counter(energy_srcs)
        metabolism_counts = Counter(metabolism)
        
        print len(orgs), 'organisms with pathway', path.name
        print 'Oxygen requirement distribution'
        print req_counts
        print 'Energy source distribution'
        print energy_counts
        print 'Metabolic category distribution'
        print metabolism_counts

if __name__ == '__main__':
    Main()
    
    
