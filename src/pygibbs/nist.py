import pylab, re, logging
from pygibbs.kegg import Kegg
from pygibbs.kegg_errors import KeggParseException
from pygibbs.thermodynamics import default_T, MissingCompoundFormationEnergy
from pygibbs.alberty import Alberty
from toolbox.util import _mkdir, calc_rmse, calc_r2
from toolbox.html_writer import HtmlWriter
from pygibbs.thermodynamic_constants import R, default_pMg
from toolbox.database import SqliteDatabase
import csv
from toolbox.ods import ODSDictReader


class NistMissingCrucialDataException(Exception):
    pass


class NistReactionBalanceException(Exception):
    pass


class NistRowData:
    def __init__(self):
        pass
    
    def ReadFromDict(self, row_dict):
        reaction = row_dict['kegg_reaction']
        if not reaction:
            raise NistMissingCrucialDataException(
                "cannot use this NIST reaction because it couldn't be mapped to KEGG IDs")

        try:
            self.sparse = NistRowData.ParseReactionFormula(reaction)
        except KeggParseException as e:
            raise NistMissingCrucialDataException("cannot use reaction \"%s\", because: %s" % (reaction, str(e)))

        if not (row_dict['K'] and row_dict['pH'] and row_dict['I']): # missing Keq, pH or I makes the data unusable
            raise NistMissingCrucialDataException(
                "cannot use this NIST reaction because it is missing information about K, pH or I")
        
        self.K = float(row_dict['K'])
        self.pH = float(row_dict['pH'])
        self.I = float(row_dict['I'])
        self.T = NistRowData.none_float(row_dict['T']) or default_T
        self.pMg = NistRowData.none_float(row_dict['pMg']) or default_pMg
        self.K_type = row_dict['Ktype']
        self.dG0_r = -R*self.T*pylab.log(self.K)
        self.evaluation = row_dict['evaluation']
        self.comment = row_dict['comment']
        self.origin = row_dict['origin']
        self.url = row_dict.get('url', None)

        if self.comment not in [None,
                                "",
                                "cosolvent => none", 
                                "added solute => none", 
                                "solvent => none"]:
            raise NistMissingCrucialDataException(
                "cannot use this NIST reaction because of an unknown solvent: " + self.comment)
            
        if self.K_type not in ["K'", "Kc'", "Km'", u"K&apos;", u"Km&apos;", u"Kc&apos;"]:
            raise NistMissingCrucialDataException(
                "cannot use this NIST reaction because of the Keq is not transformed: " + self.K_type)
    
    def ReadFromDatabase(self, row_dict):
        self.K = row_dict['K']
        self.pH = row_dict['pH']
        self.I = row_dict['I']
        self.T = row_dict['T']
        self.pMg = row_dict['pMg']
        self.K_type = row_dict['Ktype']
        self.dG0_r = row_dict['dG0']
        self.evaluation = row_dict['evaluation']
        self.comment = row_dict['comment']
        self.sparse = NistRowData.SetReactionFromString(row_dict['reaction'])
        self.origin = row_dict['origin']
        self.url = row_dict['url']
    
    @staticmethod
    def none_float(x):
        if not x:
            return None
        try:
            return float(x)
        except ValueError:
            return None
    
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
        for (cid, amount) in NistRowData.ParseReactionFormulaSide(left).iteritems():
            sparse_reaction[cid] = -amount
        for (cid, amount) in NistRowData.ParseReactionFormulaSide(right).iteritems():
            if (cid in sparse_reaction):
                raise KeggParseException("C%05d appears on both sides of this formula" % cid)
            sparse_reaction[cid] = amount
        
        return sparse_reaction
    
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
    def GetReactionString(sparse):
        return ' + '.join(["%g C%05d" % (coeff, cid) for (cid, coeff) in sorted(sparse.iteritems())])

    @staticmethod
    def SetReactionFromString(s):
        sparse = {}
        for coeff_cid in s.split(' + '):
            coeff, cid = coeff_cid.split(' ')
            coeff = float(coeff)
            cid = int(cid[1:])
            sparse[cid] = coeff
        return sparse
      
    def GetAllCids(self):
        return self.sparse.keys()
    
    def PredictFormationEnergy(self, thermodynamics, cid):
        return thermodynamics.cid2PseudoisomerMap(cid).Transform(pH=self.pH, pMg=self.pMg, I=self.I, T=self.T)
    
    def PredictReactionEnergy(self, thermodynamics):
        return sum([self.PredictFormationEnergy(thermodynamics, cid)*coeff 
                    for cid, coeff in self.sparse.iteritems()])
            
    
class Nist(object):
    def __init__(self, db, html_writer, kegg=None):
        self.db = db
        self.html_writer = html_writer
        if (kegg == None):
            self.kegg = Kegg(self.db)
        else:
            self.kegg = kegg

    def FromFile(self, filename):
        """
            Reads the contents of the CSV file into the database for faster
        """
        if filename.find('.ods') != -1:
            reader = ODSDictReader
        else:
            reader = csv.DictReader
            
        ###
        # read the NIST data from the CSV file, complete the missing conditions with the default values
        # replace the textual reaction with a "sparse_reaction" vector representation of the formula
        # and calculate the delta-G0 (instead of the Keq)
        logging.info('Reading NIST reaction data from: ' + filename)
        self.data = []
        self.cid2count = {}
        row_number = 1
        for row_dict in reader(open(filename, 'r')):
            row_number += 1
            origin = '%s - row %d' % (filename, row_number)
            row_dict['origin'] = origin
            try:
                nist_row_data = NistRowData()
                nist_row_data.ReadFromDict(row_dict)
            except NistMissingCrucialDataException as e:
                logging.debug("%s - %s" % (origin, str(e)))
                continue
            try:
                self.BalanceReaction(nist_row_data.sparse)
            except NistReactionBalanceException as e:
                logging.warning("%s - %s" % (origin, str(e)))
                logging.debug(str(nist_row_data.sparse))
                continue
            
            self.data.append(nist_row_data)
            for cid in nist_row_data.GetAllCids():
                self.cid2count[cid] = self.cid2count.setdefault(cid, 0) + 1
        
    def ToDatabase(self):
        columns = ['K REAL', 'pH REAL', 'I REAL', 'pMg REAL', 'T REAL',
                 'dG0 REAL', 'Ktype TEXT', 'evaluation TEXT', 
                 'comment TEXT', 'reaction TEXT', 'origin TEXT',
                 'url TEXT']
        self.db.CreateTable('nist_data', ','.join(columns))
        for nist_row_data in self.data:
            row_list = [nist_row_data.K, nist_row_data.pH, nist_row_data.I,
                        nist_row_data.pMg, nist_row_data.T, nist_row_data.dG0_r,
                        nist_row_data.K_type, nist_row_data.evaluation,
                        nist_row_data.comment, 
                        nist_row_data.GetReactionString(nist_row_data.sparse),
                        nist_row_data.origin,
                        nist_row_data.url]
            self.db.Insert('nist_data', row_list)
            
    def FromDatabase(self):
        self.data = []
        self.cid2count = {}
        logging.info('Reading NIST reaction data from database')
        for row_dict in self.db.DictReader('nist_data'):
            nist_row_data = NistRowData()
            nist_row_data.ReadFromDatabase(row_dict)
            self.data.append(nist_row_data)
            for cid in nist_row_data.GetAllCids():
                self.cid2count[cid] = self.cid2count.setdefault(cid, 0) + 1
        
    def Load(self):
        if not self.db.DoesTableExist('nist_data'):
            self.FromFile('../data/thermodynamics/nist.ods')
            self.ToDatabase()
        else:
            self.FromDatabase()    
    
    def BalanceReaction(self, sparse_reaction):
        """
            Checks whether a reaction is balanced.
            If there is an imbalance of oxygen or hydrogen atoms, BalanceReaction
            changes the sparse_reaction by adding H2O and H+ until it is balanced.
            
            Returns:
                True - if the reaction is balanced or if it cannot be tested at all
                False - otherwise
        """
        atom_bag = {}
        try:
            for (cid, coeff) in sparse_reaction.iteritems():
                cid_atom_bag = self.kegg.cid2atom_bag(cid)
                if (cid_atom_bag == None):
                    logging.debug("C%05d has no explicit formula, cannot check if this reaction is balanced" % cid)
                    return
                cid_atom_bag['e-'] = self.kegg.cid2compound(cid).get_num_electrons()
                
                for atomicnum, count in cid_atom_bag.iteritems():
                    atom_bag[atomicnum] = atom_bag.get(atomicnum, 0) + count*coeff
                    
        except KeyError as e:
            logging.warning(str(e) + ", cannot check if this reaction is balanced")
            return
    
        if (atom_bag.get('O', 0) != 0):
            #sys.stderr.write("WARNING: Need to add H2O to balance this reaction: " + str(sparse_reaction) + "\n")
            sparse_reaction[1] = sparse_reaction.get(1, 0) - atom_bag['O'] # balance the number of oxygens by adding C00001 (water)
            atom_bag['H'] = atom_bag.get('H', 0) - 2 * atom_bag['O'] # account for the 2 hydrogens in each added water molecule
            atom_bag['e-'] = atom_bag.get('e-', 0) - 10 * atom_bag['O'] # account for the 10 electrons in each added water molecule
            atom_bag['O'] = 0
        
        if (atom_bag.get('H', 0) != 0):
            sparse_reaction[80] = sparse_reaction.get(80, 0) - atom_bag['H'] # balance the number of hydrogens by adding C00080 (H+)
            atom_bag['H'] = 0
        
        for atomtype in atom_bag.keys():
            if atom_bag[atomtype] == 0:
                del atom_bag[atomtype]

        if atom_bag:
            raise NistReactionBalanceException("Reaction cannot be balanced: " + str(atom_bag))
    
    def GetAllCids(self):
        return sorted(self.cid2count.keys())
    
    def AnalyzeStats(self):
        """
            Produces a set of plots that show some statistics about the NIST database
        """
        logging.info('Analyzing stats against NIST database.')
        
        T_list = []
        I_list = []
        pH_list = []
        pMg_list = []
        for nist_row_data in self.data:
            T_list.append(nist_row_data.T - 273.15)
            I_list.append(nist_row_data.I)
            pH_list.append(nist_row_data.pH)
            pMg_list.append(nist_row_data.pMg)
        
        fig1 = pylab.figure()
        pylab.rc('text', usetex=False)
        pylab.rc('font', family='serif', size=7)
        pylab.title("NIST database statistics")
        pylab.subplot(2,2,1)
        pylab.hist(T_list, pylab.arange(int(min(T_list)), int(max(T_list)+1), 2.5))
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

        self.html_writer.embed_matplotlib_figure(fig1)

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

        self.html_writer.embed_matplotlib_figure(fig2)
        logging.info('Done analyzing stats.')

    def verify_results(self, thermodynamics, T_range=None):
        """Calculate all the dG0_r for the reaction from NIST and compare to
           the measured data.
        
        Write results to HTML.
        
        Args:
            thermodynamics: a Thermodynamics object that provides dG estimates.
            html_writer: to write HTML.
            ignore_I: whether or not to ignore the ionic strength in NIST.
        """
        
        known_cid_set = thermodynamics.get_all_cids()
        dG0_obs_vec = []
        dG0_est_vec = []
       
        # A mapping from each evaluation method (NIST calls separates them to
        # A, B, C and D) to the results of the relevant measurements
        evaluation_map = {}
        total_list = []
        
        cid2count = {}
        for row_data in self.data:
            for cid in row_data.GetAllCids():
                cid2count[cid] = cid2count.setdefault(cid, 0) + 1
        
        for row_data in self.data:
            unknown_set = set(row_data.GetAllCids()).difference(known_cid_set)

            if unknown_set:
                logging.debug("a compound in (%s) doesn't have a dG0_f" % row_data.origin)
                continue
            
            if T_range and not (T_range[0] < row_data.T < T_range[1]):
                logging.debug("The temperature of the measurement (%.1f) is out of range: %.1f - %.1f" %\
                              (row_data.T, T_range[0], T_range[1]) )
                continue
            
            #label = row_data.evaluation
            label = row_data.K_type
            
            if label not in evaluation_map:
                evaluation_map[label] = ([], [])
            
            try:
                dG0_pred = row_data.PredictReactionEnergy(thermodynamics)
            except MissingCompoundFormationEnergy:
                logging.debug("a compound in (%s) doesn't have a dG0_f" % row_data.origin)
                continue
                
            dG0_obs_vec.append(row_data.dG0_r)
            dG0_est_vec.append(dG0_pred)
            evaluation_map[label][0].append(row_data.dG0_r)
            evaluation_map[label][1].append(dG0_pred)
            n_measurements = min([cid2count[cid] for cid in row_data.GetAllCids()])
            error = abs(row_data.dG0_r - dG0_pred)

            total_list.append([error, row_data.dG0_r, dG0_pred, 
                               row_data.sparse, row_data.pH, row_data.pMg, 
                               row_data.I, row_data.T, row_data.evaluation, 
                               n_measurements, row_data.origin, row_data.url])
        
        if not dG0_obs_vec:
            return 0, 0
        
        # plot the profile graph
        pylab.rcParams['text.usetex'] = False
        pylab.rcParams['legend.fontsize'] = 12
        pylab.rcParams['font.family'] = 'sans-serif'
        pylab.rcParams['font.size'] = 16
        pylab.rcParams['lines.linewidth'] = 2
        pylab.rcParams['lines.markersize'] = 3
        pylab.rcParams['figure.figsize'] = [8.0, 6.0]
        pylab.rcParams['figure.dpi'] = 100
        
        fig = pylab.figure()
        pylab.hold(True)
        
        colors = ['purple', 'orange', 'lightgreen', 'red', 'cyan']
        for e in sorted(evaluation_map.keys()):
            (measured, predicted) = evaluation_map[e]
            label = '%s (N = %d, RMSE = %.2f [kJ/mol])' % (e, len(measured), calc_rmse(measured, predicted))
            c = colors.pop(0)
            pylab.plot(measured, predicted, marker='.', linestyle='None', markerfacecolor=c, markeredgecolor=c, markersize=5, label=label)
        
        pylab.legend(loc='upper left')
        
        r2 = calc_r2(dG0_obs_vec, dG0_est_vec)
        rmse = calc_rmse(dG0_obs_vec, dG0_est_vec)
        pylab.title(r'N = %d, RMSE = %.1f [kJ/mol], r$^2$ = %.2f' % (len(dG0_obs_vec), rmse, r2), fontsize=14)
        pylab.xlabel(r'$\Delta_{obs} G^\circ$ [kJ/mol]', fontsize=14)
        pylab.ylabel(r'$\Delta_{est} G^\circ$ [kJ/mol]', fontsize=14)
        min_x = min(dG0_obs_vec)
        max_x = max(dG0_obs_vec)
        pylab.plot([min_x, max_x], [min_x, max_x], 'k--')
        pylab.axis([-60, 60, -60, 60])
        self.html_writer.embed_matplotlib_figure(fig, width=400, height=300)
        
        fig = pylab.figure()
        pylab.hist([(row[1] - row[2]) for row in total_list], bins=pylab.arange(-50, 50, 0.5))
        pylab.title(r'RMSE = %.1f [kJ/mol]' % rmse, fontsize=14)
        pylab.xlabel(r'$\Delta_{obs} G^\circ - \Delta_{est} G^\circ$ [kJ/mol]', fontsize=14)
        pylab.ylabel(r'no. of measurements', fontsize=14)
        self.html_writer.embed_matplotlib_figure(fig, width=400, height=300)

        fig = pylab.figure()
        pylab.plot([row[9] for row in total_list], [abs(row[1] - row[2]) for row in total_list], '.')
        pylab.title(r'Effect of no. of measurements on estimation error', fontsize=14)
        pylab.xlabel(r'minimum no. of measurements among reaction compounds', fontsize=14)
        pylab.ylabel(r'$|| \Delta_{obs} G^\circ - \Delta_{est} G^\circ ||$ [kJ/mol]', fontsize=14)
        pylab.xscale('log')
        self.html_writer.embed_matplotlib_figure(fig, width=400, height=300)

        table_headers = ["|err|", "dG'0 (obs)", "dG'0 (est)", "reaction", "pH", "pMg", "I", "T", "eval.", "min no. measure.", "url"]
        dict_list = []
        for row in sorted(total_list, reverse=True):
            dict = {}
            dict['|err|'] = '%.1f' % row[0]
            dict['dG\'0 (obs)'] = '%.1f' % row[1]
            dict['dG\'0 (est)'] = '%.1f' % row[2]
            dict['reaction'] = self.kegg.sparse_to_hypertext(row[3], show_cids=False)
            dict['pH'] = '%.1f' % row[4]
            dict['pMg'] = '%.1f' % row[5]
            dict['I'] = '%.2f' % row[6]
            dict['T'] = '%.1f' % row[7]
            dict['eval.'] = row[8]
            dict['min no. measure.'] = row[9]
            if row[11]:
                dict['url'] = '<a href="%s">link</a>' % row[11]
            else:
                dict['url'] = ''
            dict_list.append(dict)
        self.html_writer.write_table(dict_list, table_headers)
        
        return len(dG0_obs_vec), rmse
    
    def FindRowsAccordingToReaction(self, sparse):
        rows = []
        for nist_row_data in self.data:
            if nist_row_data.sparse == sparse:
                rows.append(nist_row_data)
        return rows

if __name__ == '__main__':
    #logging.getLogger('').setLevel(logging.DEBUG)
    _mkdir("../res/nist")
    db = SqliteDatabase('../res/gibbs.sqlite')    
    html_writer = HtmlWriter("../res/nist/statistics.html")
    nist = Nist(db, html_writer)
    nist.FromFile('../data/thermodynamics/nist.ods')
    nist.ToDatabase()
    nist.AnalyzeStats()
    html_writer.close()