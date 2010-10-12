import re

STARTS_WITH_SPACE = re.compile(r'^\s')
SPLITTER = re.compile('\s+')

f = open('compounds.txt')

while True:
    line = f.readline()
    fields = SPLITTER.split(line.strip(), 1)
    if fields[0] == 'ENTRY':
        