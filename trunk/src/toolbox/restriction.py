from Bio import Entrez, SeqIO
import sys

Entrez.email = 'elad.noor@weizmann.ac.il'
handle = Entrez.efetch(db="nucleotide", rettype="fasta", id="48994873")
seq_record = SeqIO.read(handle, "fasta")
handle.close()
print "%s with %d nucleotides" % (seq_record.id, len(seq_record.seq))

restriction_enzymes = {'EcoRI':['GAATTC'],\
                       'EcoRII':['CCAGG','CCTGG'],\
                       'BamHI':['GGATCC'],\
                       'HindIII':['AAGCTT'],\
                       'TaqI':['TCGA'],\
                       'NotI':['GCGGCCGC'],\
                       'HinfI':['GAATCA','GAGTCA','GACTCA','GATTCA'],\
                       'Sau3A':['GATC'],\
                       'PovII':['CAGCTG'],\
                       'SmaI':['CCCGGG'],\
                       'HaeII':['GGCC'],\
                       'HgaI':['GACGC'],\
                       'AluI':['AGCT'],\
                       'EcoRV':['GATATC'],\
                       'KpnI':['GGTACC'],\
                       'PstI':['CTGCAG'],\
                       'SacI':['GAGCTC'],\
                       'SaII':['GTCGAC'],\
                       'ScaI':['AGTACT'],\
                       'SpeI':['ACTAGT'],\
                       'SphI':['GCATGC'],\
                       'StuI':['AGGCCT'],\
                       'SbaI':['TCTAGA']}

sys.stdout.write("%8s | %8s | %5s | %5s\n" % ("NAME", "RS", "SITES", "Average Length"))
for (name, recognition_sites) in restriction_enzymes.iteritems():
    for rs in recognition_sites:
        num_fragments = seq_record.seq.count(rs)
        sys.stdout.write("%8s | %8s | %5d | %5d\n" % (name, rs, num_fragments, len(seq_record.seq)/num_fragments))