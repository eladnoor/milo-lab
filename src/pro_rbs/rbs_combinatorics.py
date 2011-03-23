#!/usr/bin/python

import itertools
import pprint
import pylab

from rbs_calc.RBS_Calculator import RBS_Calculator

from Bio.Seq import Seq
from Bio.Alphabet import DNAAlphabet
from Bio.Align import MultipleSeqAlignment

R = 8.31e-3 # kJ/(K*mol)
JOULES_PER_CAL = 4.184
TEMP = 298.15 # K

import os
# NOTE(flamholz): set this to the location of your NUPACK install.
os.environ['NUPACKHOME'] = '/home/flamholz/Downloads/nupack3.0/'


def PredictedRate(dG):
	dG_kj = JOULES_PER_CAL * dG
	rt = R*TEMP
	return pylab.exp(-dG/rt)

def LogRate(dG):
	return pylab.log10(PredictedRate(dG))

def ExpandSeq(seq, pos):
	myseq = list(seq.tostring())
	for base in 'ATCG':
		myseq[pos] = base
		yield ''.join(myseq)

def MultiExpandSeq(seq, positions):
	myseq = list(seq.tostring())
	for bases in itertools.combinations('ATCG', len(positions)):
		for i, pos in enumerate(positions):
			myseq[pos] = bases[i]
			yield ''.join(myseq)

seq = Seq(''.join(['TCAGCAGGACGCACTGACC',
		   'GAATTCTACTAGT',
		   #'TAATAGAAATAATTTTGTTTAACTTTA',
		   #'CAACAGAAACAACCCCGCCCAATCCCA',
  		   'ACACACACACACACACACACACACACA',
		   'AGGGGATTAATTATGCATCATCACCATCACCACG']),
	  DNAAlphabet())

seq_str = seq.tostring()
print seq_str

for i, c in enumerate(seq_str):
	print '%d, %s' % (i, c)

#rbs_start = seq_str.find('AGGAGG')
rbs_start = seq_str.find('AGGGGATTAA')

#spacer_loc = seq_str.find('ACACACACATGCAT')
spacer_loc = seq_str.find('TAATTATGCATCAT')

test_positions = range(rbs_start, rbs_start + len(seq[rbs_start:spacer_loc]))
print 'Testing positions', test_positions
position_dg_std = {}
for i in test_positions:
	dgs = []
	for seq_option in ExpandSeq(seq,i):
		start_range = [0, len(seq_option)]
		name = seq_option
		calcObj = RBS_Calculator(seq_option.upper(), start_range, name)
		calcObj.calc_dG()
		if calcObj.dG_total_list:
			dgs.append(calcObj.dG_total_list[0])
		else:
			print 'Base', i
			calcObj.print_dG()

	if dgs:
		position_dg_std[i] = pylab.std(dgs)

pp = pprint.PrettyPrinter(indent=4)
print 'dG std dev by position'
pp.pprint(position_dg_std)

top_pos = sorted(position_dg_std.iteritems(), reverse=True, key=lambda x: x[1])
positions = [pos for pos, _ in top_pos[:6]]
#positions.append(top_pos[-1][0])
dg_by_seq = {}
print 'Combinatorial prediction for positions:', positions
for length in xrange(len(positions)):
	for pos_combination in itertools.combinations(positions, length):
		for seq_option in MultiExpandSeq(seq, pos_combination):
			start_range = [0, len(seq_option)]
			name = seq_option
			calcObj = RBS_Calculator(seq_option.upper(), start_range, name)
			calcObj.calc_dG()
			if calcObj.dG_total_list:
				dg_by_seq[seq_option] = calcObj.dG_total_list[0]

binding_dGs = dg_by_seq.values()
max_binding = max(binding_dGs)
min_binding = min(binding_dGs)
binding_range = max_binding - min_binding
print 'min binding dG', min_binding, 'kcal/mol, predicted rate', PredictedRate(min_binding)
print 'max binding dG', max_binding, 'kcal/mol, predicted rate', PredictedRate(max_binding)
print 'Binding range', binding_range, 'kcal/mol'

seq_by_dg = {}
for seq, dg in dg_by_seq.iteritems():
	rounded = round(dg, 1)
	seq_by_dg.setdefault(rounded, []).append(seq)
	

log_rates = [LogRate(dg) for dg in binding_dGs]
pylab.hist(log_rates)
pylab.xlabel('log10(K_binding)')
pylab.ylabel('Number of RBS')
pylab.show()

print 'sequence by dG'
#pp.pprint(seq_by_dg)

for dG in sorted(seq_by_dg.keys()):
	print '%.2g:' % dG,
	for seq in seq_by_dg[dG]:
		print seq


#print 'dG by sequence'
#pp.pprint(dg_by_seq)



