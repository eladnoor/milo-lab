import pylab, re, logging, sys
from kegg import KeggParseException, Kegg
from thermodynamics import default_T
from matplotlib.backends.backend_pdf import PdfPages
from alberty import Alberty
from toolbox import util
from toolbox.html_writer import HtmlWriter
from pygibbs.thermodynamics import default_pMg
from numpy.core.numeric import arange

class Nist:
    def __init__(self, kegg=None, fname='../data/thermodynamics/nist.csv'):
        if (kegg == None):
            self.kegg = Kegg()
        else:
            self.kegg = kegg
        self.ParseNistCsv(fname)
        
        self.cid2count = {}
        for row in self.data:
            for cid in row['sparse'].keys():
                self.cid2count[cid] = self.cid2count.setdefault(cid, 0) + 1

    def ParseNistCsv(self, filename):
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
        # read the NIST data from the CSV file, complete the missing conditions with the default values
        # replace the textual reaction with a "sparse_reaction" vector representation of the formula
        # and calculate the delta-G0 (instead of the Keq)
        self.data = []
        for row_dict in util.ReadCsvWithTitles(filename):
            if (row_dict['comment'] and row_dict['comment'] != "cosolvent => none" and row_dict['comment'] != "addded solute => none" and row_dict['comment'] != "solvent => none"):
                continue

            try:
                reaction = row_dict['kegg_reaction']
                if not reaction: # this reaction couldn't be mapped to KEGG IDs
                    continue
                row_dict['sparse'] = Nist.ParseReactionFormula(reaction)
                if not self.BalanceReaction(row_dict['sparse']):
                    continue
            except KeggParseException as e:
                logging.warning("cannot use reaction \"%s\", because: %s" % (reaction, str(e)))
                continue
    
            if not (row_dict['K'] and row_dict['pH'] and row_dict['I']): # missing Keq, pH or I makes the data unusable
                continue
            row_dict['K'] = float(row_dict['K'])
            row_dict['pH'] = float(row_dict['pH'])
            row_dict['I'] = float(row_dict['I'])
            row_dict['T'] = none_float(row_dict['T']) or default_T
            row_dict['pMg'] = none_float(row_dict['pMg']) or default_pMg
            self.data.append(row_dict)
        
    @staticmethod
    def ParseReactionFormulaSide(s):
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
    def ParseReactionFormula(formula):
        """ parse a two-sided formula such as: 2 C00001 => C00002 + C00003 
            return the set of substrates, products and the direction of the reaction
        """
        try:
            (left, right) = formula.split(' = ', 1)
        except ValueError:
            raise KeggParseException("There should be exactly one '=' sign")
        sparse_reaction = {}
        for (cid, amount) in Nist.ParseReactionFormulaSide(left).iteritems():
            sparse_reaction[cid] = -amount
        for (cid, amount) in Nist.ParseReactionFormulaSide(right).iteritems():
            if (cid in sparse_reaction):
                raise KeggParseException("C%05d appears on both sides of this formula" % cid)
            sparse_reaction[cid] = amount
        
        return sparse_reaction
    
    def BalanceReaction(self, sparse_reaction):
        atom_bag = {}
        try:
            for (cid, coeff) in sparse_reaction.iteritems():
                cid_atom_bag = self.kegg.cid2atom_bag(cid)
                if (cid_atom_bag == None):
                    logging.debug("C%05d has no explicit formula, cannot check if this reaction is balanced" % cid)
                    return True
                for (atomicnum, count) in cid_atom_bag.iteritems():
                    atom_bag[atomicnum] = atom_bag.get(atomicnum, 0) + count*coeff
                    
        except KeyError as e:
            logging.warning(str(e) + ", cannot check if this reaction is balanced")
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
    
    def AnalyzeStats(self, html_writer):
        """
            Produces a set of plots that show some statistics about the NIST database
        """

        T_list = []
        I_list = []
        pH_list = []
        pMg_list = []
        for row in self.data:
            T_list.append(float(row['T']) - 273.15)
            I_list.append(float(row['I']))
            pH_list.append(float(row['pH']))
            pMg_list.append(float(row['pMg']))
        
        fig1 = pylab.figure()
        pylab.rc('text', usetex=False)
        pylab.rc('font', family='serif', size=7)
        pylab.title("NIST database statistics")
        pylab.subplot(2,2,1)
        pylab.hist(T_list, arange(int(min(T_list)), int(max(T_list)+1), 2.5))
        pylab.xlabel("Temperature (C)")
        pylab.ylabel("No. of measurements")
        pylab.subplot(2,2,2)
        pylab.hist(pMg_list, pylab.arange(0, 6, 0.1))
        pylab.xlabel("pMg")
        pylab.ylabel("No. of measurements")
        pylab.subplot(2,2,4)
        pylab.hist(pH_list, pylab.arange(4, 11, 0.1))
        pylab.xlabel("pH")
        pylab.ylabel("No. of measurements")
        pylab.subplot(2,2,3)
        pylab.hist(I_list, pylab.arange(0, 1, 0.025))
        pylab.xlabel("Ionic Strength [mM]")
        pylab.ylabel("No. of measurements")

        html_writer.embed_matplotlib_figure(fig1)

        alberty = Alberty()
        alberty_cids = set(alberty.cid2pmap_dict.keys())
        
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
        pylab.rc('font', size=10)
        pylab.hold(True)
        p1 = pylab.bar(range(N), hist_a, color='b')
        p2 = pylab.bar(range(N), hist_b, color='r', bottom=hist_a[0:N])
        pylab.text(N-1, hist_a[N-1]+hist_b[N-1], '> %d' % (N-1), fontsize=10, horizontalalignment='right', verticalalignment='baseline')
        pylab.xlabel("N reactions")
        pylab.ylabel("no. of compounds measured in N reactions")
        pylab.legend((p1[0], p2[0]), ("Exist in Alberty's database", "New compounds"))

        html_writer.embed_matplotlib_figure(fig2)
        
if (__name__ == "__main__"):
    logging.basicConfig(level=logging.INFO, stream=sys.stderr)
    util._mkdir("../res/nist")
    nist = Nist()
    html_writer = HtmlWriter("../res/nist/statistics.html")
    nist.AnalyzeStats(html_writer)
    html_writer.close()