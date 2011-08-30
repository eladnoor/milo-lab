import itertools
import numpy as np

NucleotideMap = {'A':'T', 'C':'G', 'G':'C', 'T':'A'}

class DNAssmebler(object):

    def __init__(self, length=4):
        self.binding_mat = {}
        self.all_sense_sequences = [''.join(x) for x in 
                               itertools.product(['A','C','G','T'], repeat=length)
                               if x[-1] not in ['G','T']]
        for seq1 in self.all_sense_sequences:
            for seq2 in self.all_sense_sequences:
                self.binding_mat[seq1, seq2] = DNAssmebler._GetPairBindingEnergy(seq1, seq2)
    
    @staticmethod
    def Complement(seq):
        return ''.join([NucleotideMap[x] for x in seq])
    
    @staticmethod
    def _GetDirectBindingEnergy(seq1, seq2):
        energy = 0
        for i in xrange(len(seq1)):
            if seq1[i] == NucleotideMap[seq2[i]]:
                energy += 1
        return energy
    
    @staticmethod
    def _GetPairBindingEnergy(seq1, seq2):
        comp1 = DNAssmebler.Complement(seq1)
        comp2 = DNAssmebler.Complement(seq2)
        return max([DNAssmebler._GetDirectBindingEnergy(seq1,  seq2),
                    DNAssmebler._GetDirectBindingEnergy(seq1,  comp2),
                    DNAssmebler._GetDirectBindingEnergy(comp1, seq2),
                    DNAssmebler._GetDirectBindingEnergy(comp1, comp2)])
    
    def GetBindingEnergy(self, seq1, seq2):
        if seq1 not in self.all_sense_sequences:
            seq1 = DNAssmebler.Complement(seq1)
        if seq2 not in self.all_sense_sequences:
            seq2 = DNAssmebler.Complement(seq2)
        return self.binding_mat[seq1, seq2]
    
    def GetSubsetEnergy(self, subset):
        total_energy = 0
        for seq1, seq2 in itertools.combinations(subset, r=2):
            total_energy = max(total_energy, self.GetBindingEnergy(seq1, seq2))
        return total_energy
    
    def FindSubsetWithMinEnergy(self, size=3, fullset=None):
        best = (np.inf, None)
        if fullset is None:
            fullset = self.all_sense_sequences
        for subset in itertools.combinations(fullset, r=size):
            energy = self.GetSubsetEnergy(subset)
            curr_pair = (energy, subset)
            if energy == 2:
                return curr_pair
            best = min(best, curr_pair)
        return best
    
class Hamming(object):
    
    G4 = np.array([[1,1,1,0,0,0,0,1],
                   [1,0,0,1,1,0,0,1],
                   [0,1,0,1,0,1,0,1],
                   [1,1,0,1,0,0,1,0]])
    
    bin2nucleotide = {(0,0):'T', (0,1):'A', (1,1):'C', (1,0):'G'}
    
    bit_permutation = {0:2, 1:0, 2:4, 3:1, 4:5, 5:3, 6:6, 7:7}
    
    @staticmethod
    def Binary2Sequence(v):
        s = ""
        for i in xrange(0, len(v), 2):
            s += Hamming.bin2nucleotide[v[i],
                                        v[i+1]]
        return s
    
    def __init__(self):
        self.code = []
        for p in itertools.product([0,1], repeat=4):
            v = np.dot(Hamming.G4.T, np.array(p)) % 2
            v_perm = [v[Hamming.bit_permutation[i]] for i in xrange(len(v))]
            self.code.append(Hamming.Binary2Sequence(v_perm))
        

if __name__ == "__main__":
    assembler = DNAssmebler()
    #print assembler.FindSubsetWithMinEnergy(7)
    
    h = Hamming()
    bindmat = np.zeros((len(h.code), len(h.code)))
    for i in xrange(len(h.code)):
        for j in xrange(len(h.code)):
            bindmat[i,j] = assembler.GetBindingEnergy(h.code[i], h.code[j])
    for i in xrange(4, 10):
        print i, assembler.FindSubsetWithMinEnergy(i, h.code)
