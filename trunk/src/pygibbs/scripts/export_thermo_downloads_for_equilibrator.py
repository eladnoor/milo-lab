#!/usr/bin/python

import subprocess

thermo_sources = ('alberty', 'PGC')
phs = ['5.0', '5.5', '6.0', '6.5', '7.0', '7.5',
       '8.0', '8.5', '9.0']
i_s = '0.1'

for thermo in thermo_sources:
    for ph in phs:
        compounds_fname = 'kegg_compounds_%s_ph%s.csv' % (thermo,
                                                          ph)
        reactions_fname = 'kegg_reactions_%s_ph%s.csv' % (thermo,
                                                          ph)
        
        args = ['python',
                './pygibbs/scripts/export_kegg_data_biochemical.py',
                '-c', compounds_fname,
                '-r', reactions_fname,
                '-t', thermo,
                '-H', ph,
                '-I', i_s]
        subprocess.call(args)
