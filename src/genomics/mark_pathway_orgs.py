#!/usr/bin/python

import csv
import logging
import sys

from optparse import OptionParser

from genomics import pathway
from pygibbs.kegg import Kegg

    
def MakeOpts():
    """Returns an OptionParser object with all the default options."""
    opt_parser = OptionParser()
    opt_parser.add_option("-p", "--pathways_file",
                          dest="pathways_file",
                          help="input JSON with pathway definitions")
    opt_parser.add_option("-i", "--img_filename",
                          dest="img_filename",
                          help="input CSV file of organism data")
    opt_parser.add_option("-o", "--output_filename",
                          dest="output_filename",
                          help="output CSV file of organism data")
    opt_parser.add_option("-m", "--max_missing_genes",
                          default=0, type="int", 
                          dest="max_missing_genes",
                          help="the maximum number of pathway genes missing")
    return opt_parser


def FindOrgsWithPathway(enzyme_map, pathway_enzymes,
                        max_missing_enzymes):
    orgs_to_count = {}
    for ec_set in pathway_enzymes:
        orgs = {}
        for ec in ec_set:
            if ec not in enzyme_map:
                logging.info("Couldn't find ec %s", ec)
                continue
        
            enz = enzyme_map[ec]
            for org in enz.genes.keys():
                o = org.lower()
                orgs[o] = True
                
        for org in orgs.keys():
            orgs_to_count[org] = orgs_to_count.get(org, 0) + 1
    
    out = [k for (k,c) in orgs_to_count.iteritems()
           if len(pathway_enzymes) - c <= max_missing_enzymes]
    return sorted(out)


def ReadIMGData(filename):
    f = open(filename)
    r = csv.DictReader(f)
    rows_by_kegg = {}
    for row in r:
        kegg_id = row.get('KEGG Organism ID')
        if not kegg_id:
            continue
        
        kegg_id = kegg_id.lower()
        rows_by_kegg.setdefault(kegg_id, []).append(row)
    
    return r.fieldnames, rows_by_kegg

    

def Main():
    options, _ = MakeOpts().parse_args(sys.argv)
    assert options.pathways_file
    assert options.img_filename
    assert options.output_filename

    kegg = Kegg.getInstance(loadFromAPI=False)
    enzyme_map = kegg.ec2enzyme_map
    
    pathways = pathway.LoadPathways(options.pathways_file)
    pathway_keys = ['Has %s' % p.name for p in pathways]
    
    orgs = []
    for p in pathways:
        o = set(FindOrgsWithPathway(enzyme_map, p.enzyme_sets,
                                    options.max_missing_genes))
        print 'Found %d with %s pathway' % (len(o), p.name)
        orgs.append(o)
    
    
    fieldnames, img_data = ReadIMGData(options.img_filename)
    fieldnames.extend(pathway_keys)
    
    f = open(options.output_filename, 'w')
    w = csv.DictWriter(f, fieldnames)
    w.writeheader()
    for org, kegg_rows in img_data.iteritems():
        for pathway_key, org_set in zip(pathway_keys, orgs):
            has_pathway = org.lower() in org_set
            for row in kegg_rows:
                row[pathway_key] = has_pathway
        
        for row in kegg_rows:
            w.writerow(row)
            
    f.close()


if __name__ == '__main__':
    Main()
    print 'Done'