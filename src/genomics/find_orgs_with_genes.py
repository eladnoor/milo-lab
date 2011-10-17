#!/usr/bin/python

import csv
import sys

from optparse import OptionParser

from pygibbs.kegg import Kegg

    
def MakeOpts():
    """Returns an OptionParser object with all the default options."""
    opt_parser = OptionParser()
    opt_parser.add_option("-e", "--ecs_filename",
                          dest="ecs_filename",
                          help="input CSV file with ECs")
    opt_parser.add_option("-a", "--additional_ecs",
                          dest="additional_ecs",
                          help="comma-separated additional ECs")
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

def GetEnzymeSets(filename):
    f = open(filename)
    reader = csv.reader(f)
    sets = []
    for row in reader:
        sets.append(set(row))
    return sets


def Main():
    options, _ = MakeOpts().parse_args(sys.argv)
    assert options.ecs_filename and options.img_filename
    assert options.output_filename
    enzyme_sets = GetEnzymeSets(options.ecs_filename)
    all_ecs = reduce(lambda x,y: x.union(y), enzyme_sets)

    additional_ecs = set()
    if options.additional_ecs:
        additional_ecs.update(options.additional_ecs.split(','))
    all_ecs.update(additional_ecs)

    kegg = Kegg.getInstance(loadFromAPI=False)
    enzyme_map = kegg.ec2enzyme_map
    
    orgs_to_count = {}
    orgs_to_ec = {}
    for ec_set in enzyme_sets:
        orgs = {}
        for ec in ec_set:
            if ec not in enzyme_map:
                print 'couldnt find ec', ec
                continue
            
            enz = enzyme_map[ec]
            for org in enz.genes.keys():
                o = org.lower()
                orgs_to_ec.setdefault(o,
                    dict((e,False) for e in all_ecs))[ec] = True
                orgs[o] = True
            
        for org in orgs.keys():
            orgs_to_count[org] = orgs_to_count.get(org, 0) + 1
    
    # Count up additional ECs supplied
    for ec in additional_ecs:
        if ec not in enzyme_map:
            print 'couldnt find ec', ec
            continue
        
        enz = enzyme_map[ec]
        for org in enz.genes.keys():
            o = org.lower()
            if o not in orgs_to_ec:
                continue
            
            orgs_to_ec[o][ec] = True
    
    total = len(enzyme_sets)
    counts = [[total - count, org] for org,count in orgs_to_count.iteritems()]
    counts = sorted(filter(lambda x: x[0] <= options.max_missing_genes, counts))
    print len(counts), 'Organisms'
    print ', '.join(['%s:%d' % (o, c) for c, o in counts])
    orgs = set([org for _, org in counts])
    count_dict = dict((org, count) for count, org in counts)
    
    r = csv.DictReader(open(options.img_filename))
    missing_genes_key = 'Pathway Genes Missing'
    fieldnames = [missing_genes_key]
    fieldnames.extend(r.fieldnames)
    fieldnames.extend(list(all_ecs))
    rows_by_kegg = {}
    for row in r:
        kegg_id = row.get('KEGG Organism ID')
        if not kegg_id:
            continue
        
        kegg_id = kegg_id.lower()
        if kegg_id not in orgs:
            continue
        
        org_ecs = orgs_to_ec.get(kegg_id)
        my_row = dict(row)
        my_row.update(org_ecs)
        missing_count = count_dict.get(kegg_id)
        my_row[missing_genes_key] = missing_count
        rows_by_kegg.setdefault(kegg_id, []).append(my_row)
    
    w = csv.DictWriter(open(options.output_filename, 'w'), fieldnames)
    w.writeheader()
    for org, row_list in rows_by_kegg.iteritems():
        for row_dict in row_list:
            w.writerow(row_dict)
    
    print 'Done'
        
    
    
    

if __name__ == '__main__':
    Main()
    
    
