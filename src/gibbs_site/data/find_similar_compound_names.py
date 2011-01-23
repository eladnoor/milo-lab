#!/usr/bin/python

from util import django_utils
import csv
import logging

# NOTE(flamholz): This is crappy. We're using the real database for
# a unit test. I wish I knew of a better way.
django_utils.SetupDjango()

from gibbs import models
from gibbs import service_config

def main():
    config = service_config.Get()
    matcher = config.compound_matcher
    
    field_names = ['Compound Name', 'Compound KEGG ID',
                   'Match Name', 'Match KEGG ID', 'Match Score']
    writer = csv.DictWriter(open('name_matches.csv', 'w'),
                            field_names)
    
    seen_pairs = set()
    min_match_score = 0.7
    
    for compound in models.Compound.objects.all():
        logging.info('Finding matches for KEGG ID %s', compound.kegg_id)
        curr_id = compound.kegg_id
        for cname in compound.common_names.all():
            logging.info('Finding matches for name %s', cname.name)
            matches = matcher.Match(cname.name)
            for m in matches:
                match_id = m.value.kegg_id
                if curr_id == match_id:
                    continue
                
                pair_a = (curr_id, match_id)
                pair_b = (match_id, curr_id)
                if pair_a in seen_pairs or pair_b in seen_pairs:
                    continue
                seen_pairs.add(pair_a)
                seen_pairs.add(pair_b)
                
                if m.score < min_match_score:
                    continue
                
                row_dict = {'Compound Name': cname.name,
                            'Compound KEGG ID': curr_id,
                            'Match Name': m.key,
                            'Match KEGG ID': match_id,
                            'Match Score': m.score}
                writer.writerow(row_dict)
    

if __name__ == '__main__':
    main()
                
                
