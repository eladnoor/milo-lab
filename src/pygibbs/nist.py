import pylab, re, logging
from pygibbs.kegg import Kegg
from pygibbs.kegg_errors import KeggParseException
from pygibbs.thermodynamics import default_T, MissingCompoundFormationEnergy
from pygibbs.alberty import Alberty
from toolbox.util import _mkdir, calc_rmse, calc_r2
from toolbox.html_writer import HtmlWriter
from pygibbs.thermodynamic_constants import R, default_I, default_pMg
from toolbox.database import SqliteDatabase
import copy

class NistMissingCrucialDataException(Exception):
    pass


class NistReactionBalanceException(Exception):
    pass


class NistRowData:
    def __init__(self):
        pass
    
    @staticmethod
    def none_float(x):
        if x:
            return float(x)
        else:
            return None
    
    def ReadFromDatabase(self, row_dict):
        self.K_tag = NistRowData.none_float(row_dict['K_tag'])
        self.pH = NistRowData.none_float(row_dict['pH'])
        if not (self.K_tag and self.pH): # missing Keq or pH makes the data unusable
            raise NistMissingCrucialDataException(
                "cannot use this NIST reaction because it is missing information about K', pH or I")
        
        self.T = NistRowData.none_float(row_dict['T']) or default_T
        
        # if there is no imformation about Ionic strength of pMg, assume by
        # default that the concentration of ions is 0 (note that pMg = 14 
        # is effectively [Mg] = 0).
        self.I = NistRowData.none_float(row_dict['I']) or 0.0
        self.pMg = NistRowData.none_float(row_dict['pMg']) or 14.0 
        self.dG0_r = -R*self.T*pylab.log(self.K_tag)
        self.evaluation = row_dict['evaluation']
        self.url = row_dict['url']
        self.ref_id = row_dict['reference_id']
        reaction = row_dict['kegg_reaction']
        if not reaction:
            raise NistMissingCrucialDataException(
                "cannot use this NIST reaction because it couldn't be mapped to KEGG IDs")
        try:
            self.sparse = NistRowData.ParseReactionFormula(reaction)
        except KeggParseException as e:
            raise NistMissingCrucialDataException("cannot use reaction \"%s\", because: %s" % (reaction, str(e)))

    def GetYear(self):
        try:
            year = int(self.ref_id[0:2])
            if year < 11:
                return year + 2000
            else:
                return year + 1900
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
    def __init__(self, html_writer):
        self.db = SqliteDatabase('../data/public_data.sqlite')
        self.html_writer = html_writer
        self.kegg = Kegg.getInstance()
        self.T_range = (298, 314)
        self.override_I = None
        self.override_pMg = None
        self.FromDatabase()

    def FromDatabase(self):
        self.data = []
        self.cid2count = {}
        logging.info('Reading NIST reaction data from database')
        for row_dict in self.db.DictReader('nist_equilibrium'):
            nist_row_data = NistRowData()
            try:
                nist_row_data.ReadFromDatabase(row_dict)
                self.data.append(nist_row_data)
                for cid in nist_row_data.GetAllCids():
                    self.cid2count[cid] = self.cid2count.setdefault(cid, 0) + 1
            except NistMissingCrucialDataException:
                continue
        
    def GetAllCids(self):
        return sorted(self.cid2count.keys())
    
    def AnalyzeStats(self):
        """
            Produces a set of plots that show some statistics about the NIST database
        """
        logging.info('Analyzing stats against NIST database (%d rows).' % len(self.data))
        
        if not self.data:
            raise Exception("The database has no rows in it")
        
        T_list = []
        I_list = []
        pH_list = []
        pMg_list = []
        year_list = []
        for nist_row_data in self.data:
            pH_list.append(nist_row_data.pH)
            T_list.append(nist_row_data.T - 273.15)

            if nist_row_data.I:
                I_list.append(nist_row_data.I)
            if nist_row_data.pMg:
                pMg_list.append(nist_row_data.pMg)
            
            year = nist_row_data.GetYear()
            if year:
                year_list.append(year)
        
        fig = pylab.figure()
        pylab.title("NIST database statistics")
        pylab.hist(T_list, pylab.arange(int(min(T_list)), int(max(T_list)+1), 2.5))
        pylab.xlabel("Temperature (C)")
        pylab.ylabel("No. of measurements")
        self.html_writer.embed_matplotlib_figure(fig, width=320, height=240)

        fig = pylab.figure()
        pylab.hist(pMg_list, pylab.arange(0, 10.1, 0.1))
        pylab.xlabel("pMg")
        pylab.ylabel("No. of measurements")
        self.html_writer.embed_matplotlib_figure(fig, width=320, height=240)

        fig = pylab.figure()
        pylab.hist(pH_list, pylab.arange(4, 11, 0.1))
        pylab.xlabel("pH")
        pylab.ylabel("No. of measurements")
        self.html_writer.embed_matplotlib_figure(fig, width=320, height=240)

        fig = pylab.figure()
        pylab.hist(I_list, pylab.arange(0, 1, 0.025))
        pylab.xlabel("Ionic Strength [mM]")
        pylab.ylabel("No. of measurements")
        self.html_writer.embed_matplotlib_figure(fig, width=320, height=240)

        # histogram of publication years
        fig = pylab.figure()
        pylab.hist(year_list, pylab.arange(1930, 2010, 5))
        pylab.xlabel("Year of publication")
        pylab.ylabel("No. of measurements")
        self.html_writer.embed_matplotlib_figure(fig, width=320, height=240)

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
        
        fig = pylab.figure()
        pylab.rc('font', size=10)
        pylab.hold(True)
        p1 = pylab.bar(range(N), hist_a, color='b')
        p2 = pylab.bar(range(N), hist_b, color='r', bottom=hist_a[0:N])
        pylab.text(N-1, hist_a[N-1]+hist_b[N-1], '> %d' % (N-1), fontsize=10, horizontalalignment='right', verticalalignment='baseline')
        pylab.xlabel("N reactions")
        pylab.ylabel("no. of compounds measured in N reactions")
        pylab.legend((p1[0], p2[0]), ("Exist in Alberty's database", "New compounds"))

        self.html_writer.embed_matplotlib_figure(fig, width=320, height=240)
        logging.info('Done analyzing stats.')

    def verify_results(self, thermodynamics):
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
        
        for row_data in self.SelectRowsFromNist():
            unknown_set = set(row_data.GetAllCids()).difference(known_cid_set)
            if unknown_set:
                logging.debug("a compound in (%s) doesn't have a dG0_f" % row_data.ref_id)
                continue
            
            label = row_data.evaluation
            
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
            error = abs(row_data.dG0_r - dG0_pred)

            total_list.append([error, row_data.dG0_r, dG0_pred, 
                               row_data.sparse, row_data.pH, row_data.pMg, 
                               row_data.I, row_data.T, row_data.evaluation, 
                               row_data.url])
        
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

        table_headers = ["|err|", "dG'0 (obs)", "dG'0 (est)", "reaction", "pH", "pMg", "I", "T", "eval.", "url"]
        dict_list = []
        for row in sorted(total_list, reverse=True):
            d = {}
            d['|err|'] = '%.1f' % row[0]
            d['dG\'0 (obs)'] = '%.1f' % row[1]
            d['dG\'0 (est)'] = '%.1f' % row[2]
            d['reaction'] = self.kegg.sparse_to_hypertext(row[3], show_cids=False)
            d['pH'] = '%.1f' % row[4]
            d['pMg'] = '%.1f' % row[5]
            d['I'] = '%.2f' % row[6]
            d['T'] = '%.1f' % row[7]
            d['eval.'] = row[8]
            if row[9]:
                d['url'] = '<a href="%s">link</a>' % row[9]
            else:
                d['url'] = ''
            dict_list.append(d)
        self.html_writer.write_table(dict_list, table_headers)
        
        return len(dG0_obs_vec), rmse
    
    def SelectRowsFromNist(self, sparse=None):
        rows = []
        if sparse:
            sparse = self.kegg.BalanceReaction(sparse)
        for nist_row_data in self.data:
            if self.T_range and not (self.T_range[0] < nist_row_data.T < self.T_range[1]):
                continue 
            if sparse and nist_row_data.sparse != sparse:
                continue
            if self.override_pMg or self.override_I:
                nist_row_copy = copy.deepcopy(nist_row_data)
                if self.override_pMg:
                    nist_row_copy.pMg = self.override_pMg
                if self.override_I:
                    nist_row_copy.I = self.override_I
                rows.append(nist_row_copy)
            else:
                rows.append(nist_row_data)
        return rows   

if __name__ == '__main__':
    #logging.getLogger('').setLevel(logging.DEBUG)
    _mkdir("../res/nist")
    html_writer = HtmlWriter("../res/nist/statistics.html")
    nist = Nist(html_writer)
    nist.AnalyzeStats()
    html_writer.close()
