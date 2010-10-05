import csv

seq = csv.reader(open('../data/pro_rbs/sequences.csv', 'r'))
seq.next()
p_names = []
p_seq = []
for j in range(13):
    (id, name, sequence) = seq.next()
    p_names.append(name)
    p_seq.append(sequence)
    
r_names = []
r_seq = []
for i in range(12):
    (id, name, sequence) = seq.next()
    r_names.append(name)
    r_seq.append(sequence)    
    
seq.next()
(tmp, linker_seq, tmp2) = seq.next()
(tmp, GFP_seq, tmp2) = seq.next()
 
#result_csv = csv.writer(open('res/seq_rbs.csv', 'w'))
#for i in range(12):
#    for j in range(13):
#        name = p_names[j] + " x " + r_names[i]
#        seq = p_seq[j] + linker_seq + r_seq[i] + GFP_seq.upper()
#        pos = len(p_seq[j] + linker_seq + r_seq[i]) + 1
#        result_csv.writerow([name, seq, pos])

result_csv = csv.writer(open('../res/pro_rbs/seq_rbs.csv', 'w'))
for i in range(12):
        name = r_names[i]
        seq = linker_seq + r_seq[i] + GFP_seq.upper()
        pos = len(linker_seq + r_seq[i]) + 1
        result_csv.writerow([name, seq, pos])
