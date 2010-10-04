"""
    A short script written for Alon Wellner (alon.wellner@weizmann.ac.il).
    This generates a file containing the EC classes in KEGG and their corresponding reaction (names of substrates and products).
"""

import kegg

entry2fields_map = kegg.parse_kegg_file('../rec/reaction.txt')
ec_reaction_pairs = []
for rid in sorted(entry2fields_map.keys()):
    field_map = entry2fields_map[rid]
    if ("DEFINITION" not in field_map):
        continue
    if ("ENZYME" not in field_map):
        continue
    
    for ec in field_map["ENZYME"].split():
        reaction = "".join(field_map["DEFINITION"].split('\t'))
        ec_reaction_pairs.append(ec.split('.') + [reaction])

output = open('../res/ec_reaction_list.txt', 'w')        
for row in sorted(ec_reaction_pairs):
    output.write(".".join(row[0:4]) + '\t' + row[4] + "\n")
    
output.close()