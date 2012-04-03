import csv
import numpy as np
import matplotlib.pyplot as plt
import logging
import json

from pygibbs.thermodynamic_constants import default_T, default_pH, default_I, default_pMg,\
    symbol_df_G0, R
from pygibbs.pseudoisomer import PseudoisomerMap
from pygibbs.kegg import Kegg
from pygibbs.kegg_errors import KeggParseException,\
    KeggReactionNotBalancedException
from pygibbs.thermodynamic_errors import MissingCompoundFormationEnergy,\
    MissingReactionEnergy
from toolbox.linear_regression import LinearRegression
from toolbox.database import SqliteDatabase
import types
from pygibbs.kegg_reaction import Reaction


def GetReactionEnergiesFromFormationEnergies(S, dG0_f):
    """Calculate reaction energies from the stoichiometric matrix
       and formation energies.

    Technically, this simply performs np.dot(dG0_f, S).
    However, since some values in dG0_f might be NaN, this makes sure
    that the rows which are not affected by these NaNs are correctly
    calculated.

    Args:
        S: stoichiometric matrix - An MxN numpy.matrix
        dG0_f: formation energies - A 1xM numpy.matrix
        conc: the concentration of all compounds except H2O (default=1)
        
    Returns:
        A 1xN numpy.matrix of reaction energies (dG0_r)  
    """
    assert type(S) == np.matrix
    
    dG0_r = np.matrix(np.zeros((1, S.shape[1])))
    for r in xrange(S.shape[1]):
        nonzero_c = list(S[:, r].nonzero()[0].flat)
        dG0_r[0, r] = dG0_f[0, nonzero_c] * S[nonzero_c, r]
        
    return dG0_r

def AddConcentrationsToReactionEnergies(S, cids, T, conc):
    logc = np.ones((1, S.shape[0])) * (R * T * np.log(conc))
    if 1 in cids:
        logc[0, cids.index(1)] = 0 # H2O concentration must not change
    return logc * S


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

    def SetConditions(self, pH=None, I=None, pMg=None, T=None):
        if pH is not None:
            self.pH = pH
        if I is not None:
            self.I = I
        if pMg is not None:
            self.pMg = pMg
        if T is not None:
            self.T = T

    def SetConditionsToDefault(self):
        self.SetConditions(pH=default_pH, I=default_I, T=default_T, pMg=default_pMg)
    
    def GetConditions(self, pH=None, I=None, pMg=None, T=None):
        if pH is None:
            pH = self.pH
        if I is None:
            I = self.I
        if pMg is None:
            pMg = self.pMg
        if T is None:
            T = self.T
        return pH, I, pMg, T

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
    
    def VerifyReaction(self, reaction):
        """
            Input:
                A Reaction
            
            Raises a MissingReactionEnergy exception in case something is preventing
            this reaction from having a delta-G prediction. For example, if one of the
            compounds has a non-trivial reference point (such as guanosine=0) but that
            reference point is not balanced throughout the reaction.
        """
        missing_cids = reaction.get_cids().difference(self.get_all_cids())
        if missing_cids:
            raise MissingReactionEnergy('Some compounds have no formation energy: ' + 
                                        ', '.join(['C%05d' % cid for cid in missing_cids]),
                                        reaction.sparse)
    
    def load_bounds(self, csv_fname):
        for row in csv.DictReader(open(csv_fname, 'r')):
            cid = int(row['cid'])
            min_c = float(row['min_c'])
            max_c = float(row['max_c'])
            if cid == 0:
                self.c_range = (min_c, max_c)
            else:
                self.bounds[cid] = (min_c, max_c)
    
    def cid_to_bounds(self, cid, use_default=True):
        curr_c_range = self.bounds.get(cid, (None, None))
        if use_default:
            curr_c_range[0] = curr_c_range[0] or self.c_range[0]
            curr_c_range[1] = curr_c_range[1] or self.c_range[1]
        return curr_c_range
    
    def WriteDataToHtml(self, html_writer):
        kegg = Kegg.getInstance()
        rowdicts = []
        for cid in self.get_all_cids():
            pdata = self.cid2PseudoisomerMap(cid)
            for nH, z, nMg, dG0 in pdata.ToMatrix():
                rowdict = {}
                rowdict['KEGG ID'] = 'C%05d' % cid
                rowdict['name'] = kegg.cid2name(cid)
                rowdict['nH'] = nH
                rowdict['z'] = z
                rowdict['nMg'] = nMg
                rowdict[symbol_df_G0] = dG0
                rowdict['reference'] = pdata.GetRef(nH, z, nMg)
                if cid in self.anchors:
                    rowdict['anchor'] = 'yes'
                else:
                    rowdict['anchor'] = 'no'
                rowdicts.append(rowdict)
        
        headers = ['KEGG ID', 'name', 'nH', 'z', 'nMg', symbol_df_G0,
                   'reference', 'anchor']
        html_writer.write_table(rowdicts, headers, decimal=1)
    
    def write_data_to_csv(self, csv_fname):
        writer = csv.writer(open(csv_fname, 'w'))
        writer.writerow(['name', 'cid', 'nH', 'z', 'nMg', 'dG0'])
        for cid in sorted(self.get_all_cids()):
            name = self.kegg.cid2name(cid)
            try:
                pdata = self.cid2PseudoisomerMap(cid)
                for nH, z, nMg, dG0 in pdata.ToMatrix():
                    writer.writerow([name, "C%05d" % cid, nH, z, nMg, '%.1f' % dG0])
            except MissingCompoundFormationEnergy as e:
                logging.warning(str(e))

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
            
    def ToDatabase(self, db, table_name, error_table_name=None):
        kegg = Kegg.getInstance()
        db.CreateTable(table_name, "cid INT, nH INT, z INT, nMg INT, "
                       "dG0 REAL, compound_ref TEXT, pseudoisomer_ref TEXT, "
                       "anchor BOOL")
        if error_table_name is not None:
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
                if error_table_name is not None:
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
    
    def GetTransformedFormationEnergies(self, cids, pH=None, I=None, pMg=None, T=None):
        """ calculate the dG0_f of each compound """
        pH, I, pMg, T = self.GetConditions(pH=pH, I=I, pMg=pMg, T=T)
        
        if type(cids) == types.IntType:
            return self.cid2PseudoisomerMap(cids).Transform(pH=pH, I=I, pMg=pMg, T=T)
        elif type(cids) == types.ListType:
            dG0_f = np.matrix(np.zeros((1, len(cids))))
            for c, cid in enumerate(cids):
                try:
                    dG0_f[0, c] = self.cid2PseudoisomerMap(cid).Transform(pH=pH, I=I, pMg=pMg, T=T)
                except MissingCompoundFormationEnergy:
                    dG0_f[0, c] = np.nan
            return dG0_f
        else:
            raise ValueError("Input argument must be 'int' or 'list' of integers")
    
    def GetTransfromedKeggReactionEnergies(self, kegg_reactions,
                                           pH=None, I=None, pMg=None, T=None,
                                           conc=1):
        kegg = Kegg.getInstance()
        S, cids = kegg.reaction_list_to_S(kegg_reactions)
        return self.GetTransfromedReactionEnergies(S, cids,
                                                   pH=pH, I=I, pMg=pMg, T=T,
                                                   conc=conc)

    def GetTransfromedReactionEnergies(self, S, cids,
                                       pH=None, I=None, pMg=None, T=None,
                                       conc=1):
        """
            Returns:
                A Numpy array (column matrix) of the transformed
                formation energies of the reactions in the stiochimetric matrix
                S. The list of cids must be the same order as the columns of S.
        """
        dG0_f = self.GetTransformedFormationEnergies(cids, pH=pH, I=I, pMg=pMg, T=T)
        dG0_r = GetReactionEnergiesFromFormationEnergies(S, dG0_f)
        if conc != 1:
            dG0_r += AddConcentrationsToReactionEnergies(S, cids, T, conc)
        
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
        r2 = np.corrcoef(vec_dG0_self, vec_dG0_other)[1, 0]
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
    
    def SetConditions(self, pH=None, I=None, pMg=None, T=None):
        self.thermo[0].SetConditions(pH=pH, I=I, pMg=pMg, T=T)
        self.thermo[1].SetConditions(pH=pH, I=I, pMg=pMg, T=T)
    
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
            self.thermo[0].VerifyReaction(reaction)
            return # it's enough if the first estimator verifies this reaction
        except MissingReactionEnergy:
            pass
        
        try:
            self.thermo[1].VerifyReaction(reaction)
            return # it's enough if the second estimator verifies this reaction
        except MissingReactionEnergy:
            pass

        raise MissingReactionEnergy('None of the Thermodynamic estimators can '
                                    'calculate the Gibbs free energy of this '
                                    'reaction',
                                    reaction.sparse)

    def GetTransformedFormationEnergies(self, cids, pH=None, I=None, pMg=None, T=None):
        """
            Return the estimates of thermo[0] if all of them are known.
            Otherwise, use thermo[1] if all of them are known.
            If both have 'missing' estimates, use thermo[0] anyway.
        """
        
        dG0_f0 = self.thermo[0].GetTransformedFormationEnergies(cids, pH=pH, I=I, pMg=pMg, T=T)
        if not np.any(np.isnan(dG0_f0)):
            return dG0_f0

        dG0_f1 = self.thermo[1].GetTransformedFormationEnergies(cids, pH=pH, I=I, pMg=pMg, T=T)
        if not np.any(np.isnan(dG0_f1)):
            return dG0_f1

        return dG0_f0
        
    def GetTransfromedReactionEnergies(self, S, cids, pH=None, I=None, pMg=None, T=None, conc=1):
        """
            Find the set of reaction Gibbs energies that are completely
            consistent with thermo[0], and also close to the energies provided
            by thermo[1].
            To find this solution, we project the vector of Gibbs energies
            obtained using thermo[1] onto the subspace spanned by the columns
            of the stoichiometric matrix (where some of the values are fixed
            according to thermo[0]).
        """

        dG0_r0 = self.thermo[0].GetTransfromedReactionEnergies(S, cids, pH=pH, I=I, pMg=pMg, T=T, conc=conc)
        if np.all(np.isfinite(dG0_r0)):
            return dG0_r0
        
        # if thermo[1] cannot estimate all reactions, just use thermo[0].
        dG0_r1 = self.thermo[1].GetTransfromedReactionEnergies(S, cids, pH=pH, I=I, pMg=pMg, T=T, conc=conc)
        if np.isnan(dG0_r1).any():
            return dG0_r0

        dG0_f0 = self.thermo[0].GetTransformedFormationEnergies(cids, pH=pH, I=I, pMg=pMg, T=T, conc=conc)
        
        finite_cols = list(np.where(np.isfinite(dG0_f0))[1].flat)
        nan_cols = list(np.where(np.isnan(dG0_f0))[1].flat)
        fixed_dG0_r = dG0_f0[:, finite_cols] * S[finite_cols, :]

        P_R, P_N = LinearRegression.RowProjection(S[nan_cols, :])
        dG0_r = dG0_r1 * P_R + fixed_dG0_r * P_N
        
        return dG0_r

class ReactionThermodynamics(Thermodynamics):
    
    def __init__(self, formation_thermo, name='reaction thermodynamic data'):
        """
            arguments:
                S      - a stoichiometric matrix of the reactions from the unreliable source
                cids   - the KEGG compound IDs of the columns of S
                dG0_r  - the dG0_r' of the reactions in S, according to the unreliable source
        """
        Thermodynamics.__init__(self, name)
        self.formations = formation_thermo
        self.kegg = Kegg.getInstance()
        self.reactions = []
        self.dG0_r_primes = []
        self.cid2dG0_f = {}
        self.var_cids = []
        self.var_nullspace = None
    
    @staticmethod
    def FromCsv(csv_fname, formation_thermo):
        data = []
        for row in csv.DictReader(open(csv_fname, 'r')):
            r = Reaction.FromFormula(row['formula'])
            r.Balance(balance_water=False)
            r.SetNames(row['enzyme'])
            dG0_r_prime = float(row['dG0_r_prime'])
            pH = float(row['pH'])
            I = float(row['I'])
            pMg = float(row.get('pMg', '10'))
            T = float(row['T'])
            data.append((r, dG0_r_prime, pH, I, pMg, T))
        
        reacthermo = ReactionThermodynamics(formation_thermo)
        reacthermo.SetConditions(pH=pH, I=I, pMg=pMg, T=T)
        for r, dG0, pH, I, pMg, T in data:
            reacthermo.AddReaction(r, dG0, pH=pH, I=I, pMg=pMg, T=T)
        reacthermo._Recalculate()
        return reacthermo
    
    def AddReaction(self, kegg_reaction, dG0_r_prime, pH, I, pMg, T):
        if self.pH != pH or self.I != I or self.T != T or self.pMg != pMg:
            raise ValueError('Reverse Legendre Transform not implemented yet. '
                             'All reaction conditions must be the same')
        self.reactions.append(kegg_reaction)
        self.dG0_r_primes.append(dG0_r_prime)
         
    def AddPseudoisomer(self, cid, nH, z, nMg, dG0, ref=""):
        self.formations.AddPseudoisomer(cid, nH, z, nMg, dG0, ref)
    
    def _Recalculate(self):
        S, cids = self.kegg.reaction_list_to_S(self.reactions)
        known_cids = self.formations.get_all_cids()
        fix_rows = [i for i, cid in enumerate(cids) if cid in known_cids]
        var_rows = [i for i, cid in enumerate(cids) if cid not in known_cids]
        fix_cids = [cids[i] for i in fix_rows]
        var_cids = [cids[i] for i in var_rows]
        
        # subtract the part of the dG0 which is fixed, and leave only the part
        # which is attributed to the NaN compounds.
        fix_S = S[fix_rows, :]
        var_S = S[var_rows, :]
        
        var_P_C, var_P_L = LinearRegression.ColumnProjection(var_S)
        var_P_R, var_P_N = LinearRegression.RowProjection(var_S)

        # take all the known dG0_primes from self.formations
        dG0_f_prime = self.formations.GetTransformedFormationEnergies(fix_cids, 
                                pH=self.pH, I=self.I, T=self.T, pMg=self.pMg)

        # project the dG0_r on the column-space of var_S to eliminate inconsistencies
        # between the dG0_r and the fixed formation energies.
        # then subtract the fixed part of the dG0_r.
        var_dG0_r_prime = np.matrix(self.dG0_r_primes) * var_P_C - dG0_f_prime * fix_S
        var_dG0_f_prime, _ = LinearRegression.LeastSquares(var_S, var_dG0_r_prime)
        
        return var_cids, var_dG0_f_prime, var_P_N

    def get_all_cids(self):
        return sorted(self.formations.get_all_cids() + self.var_cids)
    
    def GetTransformedFormationEnergies(self, cids, pH=None, I=None, pMg=None, T=None):
        S = np.matrix(np.eye(len(cids)))
        return self.GetTransfromedReactionEnergies(S, cids, pH=pH, I=I, pMg=pMg, T=T)

    def GetTransfromedReactionEnergies(self, S, cids, pH=None, I=None, pMg=None, T=None, conc=1):
        if pH != None or I != None or T != None or pMg != None:
            raise MissingReactionEnergy('Cannot adjust the reaction conditions in ReactionThermodynamics', None)

        # take all the known dG0_primes from self.formations
        dG0_f_prime = self.formations.GetTransformedFormationEnergies(cids, 
                                pH=self.pH, I=self.I, T=self.T, pMg=self.pMg)

        # adjust dG0_r_primes, then perform linear regression to find the unknown dG0_f_primes
        var_cids, var_dG0_f_prime, var_P_N = self._Recalculate()

        # place the new values into the dG0_f vector at the right places (all should be NaN)
        for i, cid in enumerate(cids):
            if cid in var_cids:
                dG0_f_prime[0, i] = var_dG0_f_prime[0, var_cids.index(cid)]

        dG0_r_prime = GetReactionEnergiesFromFormationEnergies(S, dG0_f_prime)
        if conc != 1:
            dG0_r_prime += AddConcentrationsToReactionEnergies(S, cids, T, conc)

        # each row that is not orthogonal to the null-space is changed to NaN
        for j in xrange(S.shape[1]):
            v = np.matrix(np.zeros((1, len(var_cids))))
            for i, cid in enumerate(var_cids):
                if cid in cids:
                    v[0, i] = S[cids.index(cid), j]
            if (abs(v * var_P_N) > 1e-10).any():
                dG0_r_prime[0, j] = np.nan
        
        return dG0_r_prime
        
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
    
    S = np.matrix([[-1,  1,  0,  0,  0,  0,  0,  0,  0], 
                   [ 0, -1, -1,  1,  0,  0, -1,  1,  1], 
                   [ 0,  0,  0, -1,  1,  1,  0,  0,  0], 
                   [ 0, -1, -1,  0,  1,  1, -1,  1,  1]]).T
    cids = [311, 158, 10, 566, 24, 36, 2, 8, 9]
    
    print alberty.GetTransformedFormationEnergies(cids)
    print alberty.GetTransfromedReactionEnergies(S, cids)
    print pgc.GetTransfromedReactionEnergies(S, cids)
    dG0_r_primes = merged.GetTransfromedReactionEnergies(S, cids)
    print dG0_r_primes
    
