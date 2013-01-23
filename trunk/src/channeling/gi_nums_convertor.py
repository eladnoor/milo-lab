##################################################################
#
# Useful functions to convert between several representations
# of enzymes, such as gi and kegg id.
#
##################################################################

import urllib
import re

def bnum2gi (b_num):
    "Return gi from given b-number (b####) or empty string if not found"
    data = urllib.urlopen('http://rest.kegg.jp/get/eco:%s'%b_num)
    eco_re = re.search(r'NCBI-GI: (\d+)',data.read())
    data.close()
    if eco_re:
        return "gi|"+eco_re.group(1)
    return ''


def kegg2gi (kegg_id):
    "Return gi from given kegg id (k#####) or empty string if not found"
    data = urllib.urlopen('http://rest.kegg.jp/get/ko:%s'%kegg_id)
    eco_re = re.search(r'ECO: (\w+)',data.read())
    data.close()
    if eco_re:
        eco_name = eco_re.group(1)
        return bnum2gi(eco_name)
    return ''


if __name__ == '__main__':
    INPUT_FILE = 'examples/full_channeling_tabel_test.csv'
    BNUM_INDICES = (0,1) # The indices of the columns of the b-numbers
    OUTPUT_FILE = 'examples/gi_full_channeling_table_test.csv'
    
    stream = open(INPUT_FILE)
    data = map(lambda line: line.split(','), stream.readlines())
    stream.close()
    giList = map(lambda a: [bnum2gi(a[i][4:]) for i in BNUM_INDICES], data)
    stream = open(OUTPUT_FILE, 'w')
    stream.write('\n'.join(map(lambda a: ','.join(a), giList)))
    stream.close()
    print "Output was written to %s"%OUTPUT_FILE
