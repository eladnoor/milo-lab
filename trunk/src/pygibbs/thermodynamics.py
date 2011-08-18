import csv
import pylab
import logging
import json

from thermodynamic_constants import default_T, default_pH, default_I, default_pMg
from pygibbs.pseudoisomer import PseudoisomerMap
from pygibbs.kegg import Kegg
from pygibbs.kegg_errors import KeggParseException,\
    KeggReactionNotBalancedException
from pygibbs.thermodynamic_errors import MissingCompoundFormationEnergy,\
    MissingReactionEnergy
from toolbox.util import calc_r2

class Thermodynamics(object):
    def __init__(self, name="Unknown Thermodynamics"):
        self.name = name
        self.pH = default_pH
        self.I = default_I
        self.T = default_T
        self.pMg = default_pMg
        
        self.c_mid = 1e-3
        self.c_range = (1e-6, 1e-2)
        self.bounds = {}
        self.cid2source_string = {}
        self.anchors = set()

    def cid2SourceString(self, cid):
        return self.cid2source_string.get(cid, "")

    def SetPseudoisomerMap(self, cid, pmap):
        """ Add a whole Pseudoisomer Map for a compound according to its CID """
        raise NotImplementedError

    def AddPseudoisomer(self, cid, nH, z, nMg, dG0, ref=""):
        """ Add a single pseudoisomer to according to its CID """
        raise NotImplementedError
        
    def cid2PseudoisomerMap(self, cid):
        """
            Given a CID, returns the entire set of pseudoisomers (using class
            PseudoisomerMap).
        """
        raise NotImplementedError

    def get_all_cids(self):
        """ Returns a list of all CIDs that have data. """
        raise NotImplementedError

    def cid2dG0(self, cid, nH, nMg=0):
        """
            Given a CID and protonation state (number of hydrogens),
            returns the chemical dG0_f (i.e. untransformed Gibbs free energy
            of formation)
            
            One can also provide nMg (number of bound Mg2+ ions),
            but the default is 0
        """
        pmap = self.cid2PseudoisomerMap(cid)
        for p_nH, unused_p_z, p_nMg, dG0 in pmap.ToMatrix():
            if nH == p_nH and nMg == p_nMg:
                return dG0
        raise MissingCompoundFormationEnergy("C%05d doesn't have a formation "
            "energy for the species (nH=%d, nMg=%d)" % (cid, nH, nMg))

    def cid2dG0_tag(self, cid, pH=None, pMg=None, I=None, T=None):
        """
            Input:
                A CID of a compound and the aqueous conditions 
                (pH, I, pMg)
            
            Returns:
                The biochemical dG'0_f (i.e. transformed Gibbs free energy
                of formation)
        """
        pH = pH or self.pH
        I = I or self.I
        T = T or self.T
        pMg = pMg or self.pMg
        return self.cid2PseudoisomerMap(cid).Transform(pH=pH, pMg=pMg, I=I, T=T)
    
    def reaction_to_dG0(self, reaction, pH=None, pMg=None, I=None, T=None):
        """
            Input:
                A reaction in sparse representation and the aqueous conditions 
                (pH, I, pMg)
            
            Returns:
                The biochemical dG'0_r (i.e. transformed changed in Gibbs free 
                energy of reaction)
        """
        self.VerifyReaction(reaction)
        return sum([coeff * self.cid2dG0_tag(cid, pH, pMg, I, T) for 
                    cid, coeff in reaction.sparse.iteritems()])
    
    def VerifyReaction(self, reaction):
        """
            Input:
                A Reaction
            
            Raises a MissingReactionEnergy exception in case something is preventing
            this reaction from having a delta-G prediction. For example, if one of the
            compounds has a non-trivial reference point (such as guanosine=0) but that
            reference point is not balanced throughout the reaction.
        """
        missing_cids = reaction.get_cids().difference(set(self.get_all_cids()))
        if missing_cids:
            raise MissingReactionEnergy('Some compounds have no formation energy: ' + 
                                        ', '.join(['C%05d' % cid for cid in missing_cids]),
                                        reaction.sparse)
    
    def cid_to_bounds(self, cid, use_default=True):
        curr_c_min, curr_c_max = self.bounds.get(cid, (None, None))
        if not curr_c_min and use_default:
            curr_c_min = self.c_range[0]
        if not curr_c_max and use_default:
            curr_c_max = self.c_range[1]
        return (curr_c_min, curr_c_max)

    def WriteDataToHtml(self, html_writer):
        kegg = Kegg.getInstance()
        dict_list = []
        for cid in self.get_all_cids():
            for nH, z, nMg, dG0 in self.cid2PseudoisomerMap(cid).ToMatrix():
                dict = {}
                dict['cid'] = 'C%05d' % cid
                dict['name'] = kegg.cid2name(cid)
                dict['nH'] = '%d' % nH
                dict['z'] = '%d' % z
                dict['nMg'] = '%d' % nMg
                dict['dG0_f'] = '%.2f' % dG0
                if cid in self.anchors:
                    dict['anchor'] = 'yes'
                else:
                    dict['anchor'] = 'no'
                dict_list.append(dict)
        
        html_writer.write_table(dict_list, ['cid', 'name', 'nH', 'z', 
                                            'nMg', 'dG0_f', 'anchor'])
    
    def write_data_to_csv(self, csv_fname):
        writer = csv.writer(open(csv_fname, 'w'))
        writer.writerow(['cid', 'nH', 'z', 'nMg', 'dG0'])
        for cid in self.get_all_cids():
            for nH, z, nMg, dG0 in self.cid2PseudoisomerMap(cid).ToMatrix():
                writer.writerow([cid, nH, z, nMg, dG0])

    def get_json_dictionary(self):
        """Returns a JSON formatted thermodynamic data."""
        kegg = Kegg.getInstance()
        formations = []
        for cid in self.get_all_cids():
            h = {}
            h['cid'] = cid
            try:
                h['name'] = kegg.cid2name(h['cid'])
            except KeyError:
                h['name'] = None
            try:
                h['inchi'] = kegg.cid2inchi(h['cid'])
            except KeyError:
                h['inchi'] = None
            try:
                h['num_electrons'] = kegg.cid2num_electrons(h['cid'])
            except KeggParseException:
                h['num_electrons'] = None

            h['source'] = self.cid2source_string.get(cid, None)
            h['species'] = []
            for nH, z, nMg, dG0 in self.cid2PseudoisomerMap(cid).ToMatrix():
                h['species'].append({"nH":nH, "z":z, "nMg":nMg, "dG0_f":dG0})
            formations.append(h)
        
        return formations

    def write_data_to_json(self, json_fname):
        """Write JSON data to a file.
        
        Args:
            json_fname: the filename to write to.
        """
        formations = self.get_json_dictionary()
        
        json_file = open(json_fname, 'w')
        json_file.write(json.dumps(formations, indent=4))
        json_file.close()

    def override_data(self, other):
        for cid, source in other.cid2source_string.iteritems():
            self.cid2source_string[cid] = source
        
        for cid in other.get_all_cids():
            pmap = other.cid2PseudoisomerMap(cid)
            self.SetPseudoisomerMap(cid, pmap)

    def write_transformed_data_to_csv(self, csv_fname):
        writer = csv.writer(open(csv_fname, 'w'))
        writer.writerow(['cid', 'pH', 'pMg', 'I', 'T', 'dG0_tag'])
        for cid in self.get_all_cids():
            dG0_tag = self.cid2dG0_tag(cid)
            writer.writerow([cid, self.pH, self.pMg, self.I, self.T, dG0_tag])
            
    def ToDatabase(self, db, table_name, error_table_name=""):
        kegg = Kegg.getInstance()
        db.CreateTable(table_name, "cid INT, nH INT, z INT, nMg INT, "
                       "dG0 REAL, compound_ref TEXT, pseudoisomer_ref TEXT, "
                       "anchor BOOL")
        if error_table_name:
            db.CreateTable(error_table_name, 'cid INT, name TEXT, error TEXT')
        
        for cid in self.get_all_cids():
            compound_ref = self.cid2SourceString(cid)
            try:
                pmap = self.cid2PseudoisomerMap(cid)
                for nH, z, nMg, dG0 in pmap.ToMatrix():
                    pseudo_ref = pmap.GetRef(nH, z, nMg)
                    db.Insert(table_name, [cid, nH, z, nMg, dG0, compound_ref, 
                                           pseudo_ref, cid in self.anchors])
            except MissingCompoundFormationEnergy as e:
                if error_table_name:
                    db.Insert(error_table_name, [cid, kegg.cid2name(cid), str(e)])
                else:
                    logging.warning(str(e))
        db.Commit()

    def FromDatabase(self, db, table_name):
        raise NotImplementedError
        #self.cid2pmap_dict = {}
        #self.anchors = set()
        #for row in db.DictReader(table_name):
        #    self.cid2pmap_dict.setdefault(row['cid'], pseudoisomer.PseudoisomerMap())
        #    self.cid2pmap_dict[row['cid']].Add(row['nH'], row['z'], row['nMg'], 
        #                                       row['dG0'])
        #    if row['anchor']:
        #        self.anchors.add(row['cid'])
    
    def GetTransformedFormationEnergies(self, cids):
        """ calculate the dG0_f of each compound """
        dG0_f = pylab.zeros((len(cids), 1))
        for c, cid in enumerate(cids):
            try:
                dG0_f[c] = self.cid2dG0_tag(cid)
            except MissingCompoundFormationEnergy:
                # this is okay, since it means this compound's dG_f will be unbound, but only if it doesn't appear in the total reaction
                dG0_f[c] = pylab.nan
        return dG0_f
    
    def GetTransfromedReactionEnergies(self, S, cids):
        dG0_f = self.GetTransformedFormationEnergies(cids)
        dG0_r = pylab.zeros((S.shape[0], 1))
        for r in xrange(S.shape[0]):
            dG0_r[r, 0] = 0.0
            for c in pylab.find(S[r, :]):
                dG0_r[r, 0] += S[r, c] * dG0_f[c, 0]
        return dG0_r
        
    def WriteFormationEnergiesToHTML(self, html_writer, cids):
        """ calculate the dG0_f of each compound """
        kegg = Kegg.getInstance()
        html_writer.write('<table border="1">\n')
        html_writer.write('  ' + '<td>%s</td>'*6 % ("KEGG CID", "Compound Name", 
                                                    "dG0_f [kJ/mol]", "nH", "z", 
                                                    "nMg") + '\n')
        for cid in cids:
            name = kegg.cid2name(cid)
            try:
                for nH, z, nMg, dG0 in self.cid2PseudoisomerMap(cid).ToMatrix():
                    html_writer.write('<tr><td><a href="%s">C%05d</a></td><td>%s</td><td>%.2f</td><td>%d</td><td>%d</td><td>%d</td></tr>\n' % \
                                      (kegg.cid2link(cid), cid, name, dG0, nH, z, nMg))
            except MissingCompoundFormationEnergy:
                html_writer.write('<tr><td><a href="%s">C%05d</a></td><td>%s</td><td>N/A</td><td>N/A</td><td>N/A</td></tr>\n' % \
                                  (kegg.cid2link(cid), cid, name))
        html_writer.write('</table>\n')

    def CalculateCoverage(self, list_or_reactions):
        covered_counter = 0
        for reaction in list_or_reactions:
            try:
                self.VerifyReaction(reaction)
                covered_counter += 1
            except MissingReactionEnergy:
                pass
        return covered_counter
    
    def CompareOverKegg(self, html_writer, other, fig_name=None):
        """
            Compare the estimation errors of two different evaluation methods
            by calculating all the KEGG reactions which both self and other 
            can estimate, and comparing using a XY plot.
        
            Write results to HTML.
        """
        
        total_list = []
        kegg = Kegg.getInstance()
        
        for rid in sorted(kegg.get_all_rids()):
            reaction = kegg.rid2reaction(rid)
            try:
                reaction.Balance()
                dG0_self =  reaction.PredictReactionEnergy(self, 
                            pH=self.pH, pMg=self.pMg, I=self.I ,T=self.T)
                dG0_other = reaction.PredictReactionEnergy(other,
                            pH=self.pH, pMg=self.pMg, I=self.I ,T=self.T)
            except (MissingCompoundFormationEnergy, MissingReactionEnergy, 
                    KeggReactionNotBalancedException, KeyError):
                continue
                
            total_list.append({'self':dG0_self, 'other':dG0_other, 'rid':rid, 
                               'reaction':reaction})
        
        if not total_list:
            return 0, 0
        
        # plot the profile graph
        pylab.rcParams['text.usetex'] = False
        pylab.rcParams['legend.fontsize'] = 12
        pylab.rcParams['font.family'] = 'sans-serif'
        pylab.rcParams['font.size'] = 12
        pylab.rcParams['lines.linewidth'] = 2
        pylab.rcParams['lines.markersize'] = 6
        pylab.rcParams['figure.figsize'] = [6.0, 6.0]
        pylab.rcParams['figure.dpi'] = 100
        
        vec_dG0_self = pylab.array([x['self'] for x in total_list])
        vec_dG0_other = pylab.array([x['other'] for x in total_list])
        vec_rid = [x['rid'] for x in total_list]
        
        fig = pylab.figure()
        pylab.hold(True)
        max_dG0 = max(vec_dG0_self.max(), vec_dG0_other.max())
        min_dG0 = min(vec_dG0_self.min(), vec_dG0_other.min())
        pylab.plot([min_dG0, max_dG0], [min_dG0, max_dG0], 'k--', figure=fig)
        pylab.plot(vec_dG0_self, vec_dG0_other, '.', figure=fig)
        for i, rid in enumerate(vec_rid):
            pylab.text(vec_dG0_self[i], vec_dG0_other[i], '%d' % rid, fontsize=6)
        r2 = calc_r2(vec_dG0_self, vec_dG0_other)
        pylab.title("$\Delta_r G^{'\circ}$ comparison per reaction, $r^2$ = %.2f" % r2)
        pylab.xlabel(self.name + ' (in kJ/mol)', figure=fig)
        pylab.ylabel(other.name + ' (in kJ/mol)', figure=fig)
        html_writer.embed_matplotlib_figure(fig, width=200, height=200, name=fig_name)
        
class ThermodynamicsWithCompoundAbundance(Thermodynamics):
    
    def __init__(self, name="Unknown ThermodynamicsWithCompoundAbundance"):
        Thermodynamics.__init__(self, name)
        self.media_list = []
        self.default_c0 = 1 # Molar
        self.cid2conc = {}
    
    def LoadConcentrationsFromBennett(self, filename=
            '../data/thermodynamics/compound_abundance.csv'):
        
        self.media_list = ['glucose', 'glycerol', 'acetate']
        for row in csv.DictReader(open(filename, 'r')):
            if not row['KEGG ID']:
                continue
            cid = int(row['KEGG ID'])
            if row['use'] != '1':
                continue
            for media in self.media_list:
                try:
                    conc = float(row[media])
                    self.cid2conc[(cid, media)] = conc
                except ValueError:
                    pass

    def LoadConcentrationsFromDB(self, db, table_name='compound_abundance'):
        self.media_list = []
        if db.DoesTableExist(table_name):
            for row in db.Execute("SELECT media FROM compound_abundance GROUP BY media"):
                self.media_list.append(row[0])
                
            self.cid2conc = {}
            for row in self.db.Execute("SELECT cid, media, concentration FROM compound_abundance"):
                cid, media, conc = row
                self.cid2conc[(cid, media)] = conc # in [M]

    def GetConcentration(self, cid, c0=None, media=None):
        c0 = c0 or self.default_c0
        if cid == 1: # the concentration of water must always be 1 [M]
            return 1
        if not media:
            return c0
        return self.cid2conc.get((cid, media), c0)
    
class PsuedoisomerTableThermodynamics(ThermodynamicsWithCompoundAbundance):
    
    def __init__(self, name="Unknown PsuedoisomerTableThermodynamics"):
        ThermodynamicsWithCompoundAbundance.__init__(self, name)
        self.cid2pmap_dict = {}
    
    @staticmethod
    def _FromDictReader(reader, thermo, label=None, 
                        name="Unknown PsuedoisomerTableThermodynamics", 
                        warn_for_conflicting_refs=True):
        """Internal and for testing only.
        
        Creates a Thermodynamics object from a DictReader.
        
        Arguments:
            reader - a list of dictionaries containing the row data
            thermo - the Thermodynamics class into which to write the data
            label - if not None, only rows with this label will be used
                    otherwise, all rows except those with the label 'skip' will be used
            warn_for_conflicting_refs - if True, print warnings if two rows
                    with the same CID have different refs.
        """
        for row in reader:
            if label and row.get('use for', None) != label:
                continue
            elif row.get('use for', None) == 'skip':
                continue
            
            cid = int(row['cid'])
            if not cid:
                continue
            
            try:
                nH = int(row['nH'])
                z = int(row['z'])
                nMg = int(row['nMg'])
                dG0 = float(row['dG0'])
                compound_ref = row.get('compound_ref', '')
                pseudo_ref = row.get('pseudoisomer_ref', '')
                thermo.AddPseudoisomer(cid, nH, z, nMg, dG0, pseudo_ref)
                if cid in thermo.cid2source_string and thermo.cid2source_string[cid] != compound_ref:
                    if warn_for_conflicting_refs:
                        logging.warning('There are conflicting references for C%05d ' % cid)
                else:
                    thermo.cid2source_string[cid] = compound_ref
            except ValueError as e:
                raise ValueError("For C%05d: %s" % (cid, str(e)))
            
        return thermo
    
    @staticmethod
    def FromDatabase(db, table_name, label=None, name=None):
        """
            Imports the pseudoisomer maps from a CSV file, with these headers:
            'cid', 'nH', 'z', 'nMg', 'dG0'
        """
        if name is None:
            name = "From table " + table_name
        reader = db.DictReader(table_name)
        thermo = PsuedoisomerTableThermodynamics(name)
        PsuedoisomerTableThermodynamics._FromDictReader(
                                                reader, thermo, label, name)
        return thermo

    @staticmethod
    def FromCsvFile(filename, label=None, name=None):
        """
            Imports the pseudoisomer maps from a CSV file, with these headers:
            'cid', 'nH', 'nMg', 'dG0'
        """
        if name is None:
            name = "From file " + filename
        reader = csv.DictReader(open(filename, 'r'))
        thermo = PsuedoisomerTableThermodynamics(name)
        PsuedoisomerTableThermodynamics._FromDictReader(
                                                reader, thermo, label, name)
        return thermo

    def cid2PseudoisomerMap(self, cid):
        if cid in self.cid2pmap_dict:
            return self.cid2pmap_dict[cid]
        else:
            format_str = ("The compound C%05d does not have a value "
                          "for its formation energy of any of its pseudoisomers")
            raise MissingCompoundFormationEnergy(format_str % cid, cid)

    def SetPseudoisomerMap(self, cid, pmap):
        self.cid2pmap_dict[cid] = pmap

    def AddPseudoisomer(self, cid, nH, z, nMg, dG0, ref=''):
        self.cid2pmap_dict.setdefault(cid, PseudoisomerMap())
        self.cid2pmap_dict[cid].Add(nH, z, nMg, dG0, ref)

    def get_all_cids(self):
        return sorted(self.cid2pmap_dict.keys())

    def test(self):
        pMg_vec = pylab.arange(0, 10, 0.1)
        dG_vec = []
        dG_f_mat = []
        sparse = {20:-1, 13:-1, 147:1, 119:1}
        #sparse = {2:-1, 1:-1, 20:1, 13:1}
        #sparse = {2:-1, 1:-1, 8:1, 9:1}
        for pMg in pMg_vec:
            dG = 0
            dG_f_vec = []
            for cid, coeff in sparse.iteritems():
                dG_f = self.cid2PseudoisomerMap(cid).Transform(pH=7.4, 
                    pMg=pMg, I=0.0, T=303.1)
                dG_f_vec.append(dG_f)
                dG += coeff * dG_f
            dG_f_mat.append(dG_f_vec)
            dG_vec.append(dG)
        pylab.plot(pMg_vec, dG_vec)
        #pylab.plot(pMg_vec, dG_f_mat)
        pylab.show()

if __name__ == "__main__":
    T = PsuedoisomerTableThermodynamics.FromCsvFile(
        '../data/thermodynamics/alberty_pseudoisomers.csv')
    #T.test()
    all_reactions = Kegg.getInstance().AllReactions()
    print "Alberty coverage over KEGG: ", T.CalculateCoverage(all_reactions)
    print "Out of %d reactions" % len(all_reactions)