import csv, pylab, re, sys
from kegg import KeggParseException, Kegg
from thermodynamics import default_T
from matplotlib.backends.backend_pdf import PdfPages
from alberty import Alberty

class Nist:
    def __init__(self, kegg=None, fname='../data/thermodynamics/nist.csv'):
        if (kegg == None):
            self.kegg = Kegg()
        else:
            self.kegg = kegg
        self.parse_nist_csv(fname)
        
        self.cid2count = {}
        for row in self.data:
            for cid in row[6].keys():
                self.cid2count[cid] = self.cid2count.setdefault(cid, 0) + 1

    def parse_nist_csv(self, fname):
        """
            Reads the table of NIST reaction measurements and returns it as a table
        """
        ###
        def none_float(x):
            try:
                return float(x)
            except ValueError:
                return None
        ###
        
        csv_reader = csv.reader(open(fname, 'r'))
        
        # read the NIST data from the CSV file, complete the missing conditions with the default values
        # replace the textual reaction with a "sparse_reaction" vector representation of the formula
        # and calculate the delta-G0 (instead of the Keq)
        self.titles = csv_reader.next()
        self.data = []
        rowid = 0
        for row in csv_reader:
            rowid += 1
            if (row[12] != "" and row[12] != "cosolvent => none" and row[12] != "addded solute => none" and row[12] != "solvent => none"):
                continue

            try:
                reaction = row[6]
                if (reaction == ""): # this reaction couldn't be mapped to KEGG IDs
                    continue
                sparse_reaction = Nist.parse_reaction_formula(reaction)
                if (not self.balance_reaction(sparse_reaction)):
                    continue
                row[6] = sparse_reaction
            except KeggParseException as e:
                sys.stderr.write("WARNING: Cannot use reaction \"%s\", because: %s\n" % (reaction, str(e)))
                continue
    
            [Keq, T, I, pH] = [none_float(x) for x in row[8:12]]
            if (Keq == None or pH == None or I == None): # missing Keq or pH makes the data unusable
                continue
            if (T == None):
                T = default_T # default temperature (room)
            row[8] = float(Keq)
            row[9] = T
            row[10] = I
            row[11] = pH
            self.data.append(row)
        
    @staticmethod
    def parse_reaction_formula_side(s):
        """ parse the side formula, e.g. '2 C00001 + C00002 + 3 C00003'
            return the set of CIDs, ignore stoichiometry
        """
        compound_bag = {}
        for member in re.split('\s+\+\s+', s):
            tokens = member.split(None, 1)
            if (len(tokens) == 1):
                amount = 1
                cid = member
            else:
                amount = tokens[0]
                cid = tokens[1]
            if (cid[0] != 'C'):
                raise KeggParseException("Compound ID does not start with a C: " + cid)
            try:
                cid = int(cid[1:])
            except ValueError:
                raise KeggParseException("Compound ID is not an integer number: " + cid)
            compound_bag[cid] = compound_bag.get(cid, 0) + float(amount)
        
        return compound_bag
      
    @staticmethod
    def parse_reaction_formula(formula):
        """ parse a two-sided formula such as: 2 C00001 => C00002 + C00003 
            return the set of substrates, products and the direction of the reaction
        """
        try:
            (left, right) = formula.split(' = ', 1)
        except ValueError:
            raise KeggParseException("There should be exactly one '=' sign")
        sparse_reaction = {}
        for (cid, amount) in Nist.parse_reaction_formula_side(left).iteritems():
            sparse_reaction[cid] = -amount
        for (cid, amount) in Nist.parse_reaction_formula_side(right).iteritems():
            if (cid in sparse_reaction):
                raise KeggParseException("C%05d appears on both sides of this formula" % cid)
            sparse_reaction[cid] = amount
        
        return sparse_reaction
    
    def balance_reaction(self, sparse_reaction):
        atom_bag = {}
        try:
            for (cid, coeff) in sparse_reaction.iteritems():
                cid_atom_bag = self.kegg.cid2atom_bag(cid)
                if (cid_atom_bag == None):
                    sys.stderr.write("WARNING: C%05d has no explicit formula, cannot check if this reaction is balanced\n" % cid)
                    return True
                for (atomicnum, count) in cid_atom_bag.iteritems():
                    atom_bag[atomicnum] = atom_bag.get(atomicnum, 0) + count*coeff
                    
        except KeyError as e:
            sys.stderr.write("WARNING: " + str(e) + ", cannot check if this reaction is balanced\n")
            return True
    
        if (atom_bag.get('O', 0) != 0):
            #sys.stderr.write("WARNING: Need to add H2O to balance this reaction: " + str(sparse_reaction) + "\n")
            sparse_reaction[1] = sparse_reaction.get(1, 0) - atom_bag['O'] # balance the number of oxygens by adding C00001 (water)
            atom_bag[1] = atom_bag.get(1, 0) - 2 * atom_bag['O'] # account for the hydrogens in the added water molecules
            atom_bag[8] = 0
        
        if (atom_bag.get('H', 0) != 0):
            sparse_reaction[80] = sparse_reaction.get(80, 0) - atom_bag['H'] # balance the number of hydrogens by adding C00080 (H+)
            atom_bag['H'] = 0
            
        for (atomicnum, balance) in atom_bag.iteritems():
            if (balance != 0):
                return False
        
        return True
    
    def analyze_stats(self, pdf_fname="../res/nist_statistics.pdf"):
        """
            Produces a set of plots that show some statistics about the NIST database
        """

        T_list = []
        I_list = []
        pH_list = []
        for row in self.data:
            T_list.append(float(row[9]) - 273.15)
            I_list.append(float(row[10]))
            pH_list.append(float(row[11]))
        
        fig1 = pylab.figure()
        pylab.rc('text', usetex=False)
        pylab.rc('font', family='serif', size=7)
        pylab.title("NIST database statistics")
        pylab.subplot(2,2,1)
        pylab.hist(T_list, range(int(min(T_list)), int(max(T_list)+1)))
        pylab.xlabel("Temperature (C)")
        pylab.ylabel("No. of measurements")
        pylab.subplot(2,2,2)
        pylab.plot(pH_list, T_list, '.', markersize=2)
        pylab.xlabel("pH")
        pylab.ylabel(r"Temperature (C)")
        pylab.subplot(2,2,4)
        pylab.hist(pH_list, pylab.arange(4, 11, 0.1))
        pylab.xlabel("pH")
        pylab.ylabel("No. of measurements")
        pylab.subplot(2,2,3)
        pylab.hist(I_list, pylab.arange(0, 1, 0.025))
        pylab.xlabel("Ionic Strength [mM]")
        pylab.ylabel("No. of measurements")

        alberty = Alberty()
        alberty_cids = set(alberty.cid2pmap.keys())
        
        N = 60 # cutoff for the number of counts in the histogram
        hist_a = pylab.zeros(N)
        hist_b = pylab.zeros(N)
        for (cid, cnt) in self.cid2count.iteritems():
            if (cnt >= N):
                cnt = N-1
            if (cid in alberty_cids):
                hist_a[cnt] += 1
            else:
                hist_b[cnt] += 1
        hist_a[0] = len(alberty_cids.difference(self.cid2count.keys()))
        
        fig2 = pylab.figure()
        pylab.rc('text', usetex=True)
        pylab.rc('font', size=10)
        pylab.hold(True)
        p1 = pylab.bar(range(N), hist_a, color='b')
        p2 = pylab.bar(range(N), hist_b, color='r', bottom=hist_a[0:N])
        pylab.text(N-1, hist_a[N-1]+hist_b[N-1], '$\ge$%d' % (N-1), fontsize=10, horizontalalignment='right', verticalalignment='baseline')
        pylab.xlabel("$N$ reactions")
        pylab.ylabel("no. of compounds measured in $N$ reactions")
        pylab.legend((p1[0], p2[0]), ("Exist in Alberty's database", "New compounds"))

        pp = PdfPages(pdf_fname)
        pp.savefig(fig1)
        pp.savefig(fig2)
        pp.close()
        
if (__name__ == "__main__"):
    nist = Nist()
    nist.analyze_stats()
