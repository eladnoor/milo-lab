import csv
import numpy as np
import matplotlib.pyplot as plt
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
from toolbox.linear_regression import LinearRegression
from toolbox.database import SqliteDatabase

def GetReactionEnergiesFromFormationEnergies(S, dG0_f):
    """
        Technically, this simply performs np.dot(S, dG0_f).
        However, since some values in dG0_f might be NaN, this makes sure
        that the rows which are not affected by these NaNs are correctly
        calculated.
    """
    dG0_r = np.zeros((S.shape[0], 1))
    for r in xrange(S.shape[0]):
        dG0_r[r, 0] = 0.0
        for c in np.nonzero(S[r, :])[0]:
            dG0_r[r, 0] += S[r, c] * dG0_f[c, 0]
    return dG0_r    

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

    def SetConditions(self, pH=None, I=None, T=None, pMg=None):
        self.pH = pH or self.pH
        self.I = I or self.I
        self.T = T or self.T
        self.pMg = pMg or self.pMg

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

    def GetJSONDictionary(self):
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

    def WriteDataToJSON(self, json_fname):
        """Write JSON data to a file.
        
        Args:
            json_fname: the filename to write to.
        """
        formations = self.GetJSONDictionary()
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
        dG0_f = np.zeros((len(cids), 1))
        for c, cid in enumerate(cids):
            try:
                dG0_f[c] = self.cid2dG0_tag(cid)
            except MissingCompoundFormationEnergy:
                # this is okay, since it means this compound's dG_f will be unbound, but only if it doesn't appear in the total reaction
                dG0_f[c] = np.nan
        return dG0_f
    
    def GetTransfromedReactionEnergies(self, S, cids):
        """
            Returns:
                A Numpy array (column matrix) of the transformed
                formation energies of the reactions in the stiochimetric matrix
                S. The list of cids must be the same order as the columns of S.
        """
        dG0_f = self.GetTransformedFormationEnergies(cids)
        return GetReactionEnergiesFromFormationEnergies(S, dG0_f)
        
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
        plt.rcParams['text.usetex'] = False
        plt.rcParams['legend.fontsize'] = 12
        plt.rcParams['font.family'] = 'sans-serif'
        plt.rcParams['font.size'] = 12
        plt.rcParams['lines.linewidth'] = 2
        plt.rcParams['lines.markersize'] = 6
        plt.rcParams['figure.figsize'] = [6.0, 6.0]
        plt.rcParams['figure.dpi'] = 100
        
        vec_dG0_self = np.array([x['self'] for x in total_list])
        vec_dG0_other = np.array([x['other'] for x in total_list])
        vec_rid = [x['rid'] for x in total_list]
        
        fig = plt.figure()
        fig.hold(True)
        max_dG0 = max(vec_dG0_self.max(), vec_dG0_other.max())
        min_dG0 = min(vec_dG0_self.min(), vec_dG0_other.min())
        plt.plot([min_dG0, max_dG0], [min_dG0, max_dG0], 'k--', figure=fig)
        plt.plot(vec_dG0_self, vec_dG0_other, '.', figure=fig)
        for i, rid in enumerate(vec_rid):
            plt.text(vec_dG0_self[i], vec_dG0_other[i], '%d' % rid, fontsize=6)
        r2 = calc_r2(vec_dG0_self, vec_dG0_other)
        plt.title("$\Delta_r G^{'\circ}$ comparison per reaction, $r^2$ = %.2f" % r2)
        plt.xlabel(self.name + ' (in kJ/mol)', figure=fig)
        plt.ylabel(other.name + ' (in kJ/mol)', figure=fig)
        html_writer.embed_matplotlib_figure(fig, width=200, height=200, name=fig_name)
        
   
class PsuedoisomerTableThermodynamics(Thermodynamics):
    
    def __init__(self, name="Unknown PsuedoisomerTableThermodynamics"):
        Thermodynamics.__init__(self, name)
        self.cid2pmap_dict = {}

    def Clone(self):
        other = PsuedoisomerTableThermodynamics()
        other.name = self.name
        other.pH = self.pH
        other.I = self.I
        other.T = self.T
        other.pMg = self.pMg
        other.c_mid = self.c_mid
        other.c_range = self.c_range
        other.bounds = dict(self.bounds)
        other.cid2source_string = dict(self.cid2source_string)
        other.anchors = set(self.anchors)
        other.cid2pmap_dict = dict(self.cid2pmap_dict)
        return other
    
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
        pMg_vec = np.arange(0, 10, 0.1)
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
        plt.plot(pMg_vec, dG_vec)
        #pylab.plot(pMg_vec, dG_f_mat)
        plt.show()

class BinaryThermodynamics(Thermodynamics):
    
    def __init__(self, thermo0, thermo1):
        Thermodynamics.__init__(self, name=thermo0.name + ' / ' + thermo1.name)
        self.thermo = (thermo0, thermo1)

    def cid2PseudoisomerMap(self, cid):
        for thermo in self.thermo:
            try:
                return thermo.cid2PseudoisomerMap(cid)
            except MissingCompoundFormationEnergy:
                continue
        format_str = ("The compound C%05d does not have a value "
                      "for its formation energy of any of its pseudoisomers")
        raise MissingCompoundFormationEnergy(format_str % cid, cid)
    
    def AddPseudoisomer(self, cid, nH, z, nMg, dG0, ref=""):
        self.thermo[0].AddPseudoisomer(cid, nH, z, nMg, dG0, ref)
        self.thermo[1].AddPseudoisomer(cid, nH, z, nMg, dG0, ref)
    
    def get_all_cids(self):
        cids = set(self.thermo[0].get_all_cids() + self.thermo[1].get_all_cids())
        return sorted(cids)
    
    def reaction_to_dG0(self, reaction, pH=None, pMg=None, I=None, T=None):
        """
            Input:
                A reaction in sparse representation and the aqueous conditions 
                (pH, I, pMg)
            
            Returns:
                The biochemical dG'0_r (i.e. transformed changed in Gibbs free 
                energy of reaction)
        """
        for thermo in self.thermo:
            try:
                # if at least one of the 'thermodynamics' in the stack verify
                # this reaction, then it is okay.
                return thermo.reaction_to_dG0(reaction, pH, pMg, I, T)
            except MissingReactionEnergy:
                continue
        
        raise MissingReactionEnergy('None of the Thermodynamic estimators can '
                                    'calculate the Gibbs free energy of this '
                                    'reaction')
            
    def SetConditions(self, pH=None, I=None, T=None, pMg=None):
        self.thermo[0].SetConditions(pH, I, T, pMg)
        self.thermo[1].SetConditions(pH, I, T, pMg)
    
    def VerifyReaction(self, reaction):
        """
            Input:
                A Reaction
            
            Raises a MissingReactionEnergy exception in case something is preventing
            this reaction from having a delta-G prediction. For example, if one of the
            compounds has a non-trivial reference point (such as guanosine=0) but that
            reference point is not balanced throughout the reaction.
        """
        try:
            # if at least one of the 'thermodynamics' in the stack verify
            # this reaction, then it is okay.
            self.thermo[0].VerifyReaction(reaction)
            self.thermo[1].VerifyReaction(reaction)
        except MissingReactionEnergy:
            raise MissingReactionEnergy('None of the Thermodynamic estimaters can '
                                        'calculate the Gibbs free energy of this '
                                        'reaction')

    def GetTransformedFormationEnergies(self, cids):
        """
            Return the estimates of thermo[0] if all of them are known.
            Otherwise, use thermo[1] if all of them are known.
            If both have 'missing' estimates, use thermo[0] anyway.
        """
        
        dG0_f0 = self.thermo[0].GetTransformedFormationEnergies(cids)
        if not np.any(np.isnan(dG0_f0)):
            return dG0_f0

        dG0_f1 = self.thermo[1].GetTransformedFormationEnergies(cids)
        if not np.any(np.isnan(dG0_f1)):
            return dG0_f1

        return dG0_f0
    
    def GetTransfromedReactionEnergies(self, S, cids):
        """
            Find the set of reaction Gibbs energies that are completely
            consistent with thermo[0], and also close to the energies provided
            by thermo[1].
            To find this solution, we project the vector of Gibbs energies
            obtained using thermo[1] onto the subspace spanned by the columns
            of the stoichiometric matrix (where some of the values are fixed
            according to thermo[0]).
        """

        dG0_r0 = self.thermo[0].GetTransfromedReactionEnergies(S, cids)
        if np.all(np.isfinite(dG0_r0)):
            return dG0_r0
        
        # if thermo[1] cannot estimate all reactions, just use thermo[0].
        dG0_r1 = self.thermo[1].GetTransfromedReactionEnergies(S, cids)
        if np.any(np.isnan(dG0_r1)):
            return dG0_r0

        dG0_f0 = self.thermo[0].GetTransformedFormationEnergies(cids)
        
        finite_cols = np.where(np.isfinite(dG0_f0))[0]
        nan_cols = np.where(np.isnan(dG0_f0))[0]
        fixed_dG0_r = np.dot(S[:, finite_cols], dG0_f0[finite_cols])

        P_C, P_L = LinearRegression.ColumnProjection(S[:, nan_cols])
        dG0_r = np.dot(P_C, dG0_r1) + np.dot(P_L, fixed_dG0_r)
        
        return dG0_r
        
if __name__ == "__main__":

    from pygibbs.groups import GroupContribution
    
    db_public = SqliteDatabase('../data/public_data.sqlite')
    db_gibbs = SqliteDatabase('../res/gibbs.sqlite')
    alberty = PsuedoisomerTableThermodynamics.FromDatabase(\
                        db_public, 'alberty_pseudoisomers', name='alberty')
    
    pgc = GroupContribution(db=db_gibbs, transformed=False)
    pgc.init()
    pgc.name = "PGC"
    
    merged = BinaryThermodynamics(alberty, pgc)
    
    S = np.array([[-1, 1, 0, 0, 0, 0, 0, 0, 0], [0, -1, -1, 1, 0, 0, -1, 1, 1], [0, 0, 0, -1, 1, 1, 0, 0, 0], [0, -1, -1, 0, 1, 1, -1, 1, 1]])
    cids = [311, 158, 10, 566, 24, 36, 2, 8, 9]
    
    print alberty.GetTransfromedReactionEnergies(S, cids).T
    print pgc.GetTransfromedReactionEnergies(S, cids).T
    print merged.GetTransfromedReactionEnergies(S, cids).T