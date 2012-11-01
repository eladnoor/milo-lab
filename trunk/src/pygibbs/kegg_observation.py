import json, logging, csv
import numpy as np
from pygibbs.kegg import Kegg
from pygibbs.thermodynamics import PsuedoisomerTableThermodynamics
from pygibbs.thermodynamic_constants import default_pH, default_T,\
    RedoxCarriers, symbol_dr_G0_prime, symbol_dr_G0, default_pMg, R
from pygibbs.nist import Nist
from pygibbs.dissociation_constants import MissingDissociationConstantError,\
    DissociationConstants
from toolbox.html_writer import NullHtmlWriter
from pygibbs.kegg_reaction import Reaction
from pygibbs import thermodynamic_constants
from toolbox.molecule import Molecule

class KeggObservation(object):
    
    TYPE_FORMATION = 1
    TYPE_REACTION = 2
    TYPE_ACID_BASE = 3
    TYPE_MG = 4
    TYPE_REDOX = 5
    
    ATOM2ELEMENT = {'C': ('Carbon', 1),
                    'H': ('H2(g)', 0.5),
                    'O': ('O2(g)', 0.5),
                    'N': ('N2(g)', 0.5),
                    'P': ('Phosphorus', 1),
                    'S': ('Sulfur', 1),
                    'Cl': ('Cl2(g)', 0.5),
                    'F': ('F2(g)', 0.5),
                    'Fe': ('Iron', 1),
                    'Mn': ('Manganese', 1),
                    'Mg': ('Magnesium', 1),
                    'Br': ('Br2(l)', 0.5),
                    'I': ('I2(s)', 0.5)}

    
    def __init__(self, obs_id, obs_type, url, anchored, dG0, sparse):
        self.obs_id = obs_id
        self.obs_type = obs_type
        self.url = url
        self.dG0 = dG0
        self.sparse = sparse
        self.anchored = anchored

    def ToDatabase(self, db, table_name):
        db.Insert(table_name, self.ToCsvRow())
    
    @staticmethod 
    def FromDatabaseRow(row):
        sparse = {}
        for cid, coeff in json.loads(row['reaction']).iteritems():
            sparse[int(cid)] = float(coeff)
        
        return KeggObservation(row['id'], row['type'], row['url'],
                                row['anchored'], row['dG0'], sparse)
        
    def ToCsvRow(self):
        return [self.obs_id, self.obs_type, self.url, self.anchored, self.dG0, 
                json.dumps(self.sparse)]
    
    def __hash__(self):
        return json.dumps(self.sparse)
    
    def __eq__(self, other):
        return self.sparse == other.sparse
    
    def Normalize(self):
        """
            In self.sparse, the first non-zero coefficient must be negative,
            where the order is defined alphabetically by the keys.
        """
        if self.sparse[min(self.sparse.keys())] > 0:
            self.sparse = dict((k,-v) for (k,v) in self.sparse.iteritems())
            self.dG0 = -self.dG0
        
class KeggObervationCollection(object):
    
    def __init__(self, html_writer, dissociation, transformed=False,
                 pH=default_pH, I=0, pMg=14, T=default_T):
        self.pH = pH
        self.I = I
        self.pMg = pMg
        self.T = T

        self.kegg = Kegg.getInstance()
        self.transformed = transformed
        self.dissociation = dissociation
        if self.dissociation is not None:
            self.cid2nH_nMg = self.dissociation.GetCid2nH_nMg(
                                    pH=self.pH, I=self.I, pMg=self.pMg, T=self.T)
        self.observations = []
        self.html_writer = html_writer
        if transformed:
            self.gibbs_symbol = symbol_dr_G0_prime
        else:
            self.gibbs_symbol = symbol_dr_G0
            
        self.FormationEnergyFileName = "../data/thermodynamics/formation_energies.csv"

    @staticmethod
    def FromFiles(html_writer, dissociation, transformed=False,
                  pH=default_pH, I=0, pMg=14, T=default_T,
                  formation_energy_fname=None,
                  anchored_redox=True):
        
        obs_collections = KeggObervationCollection(html_writer,
            dissociation, transformed, pH=pH, I=I, pMg=pMg, T=T)
        
        if formation_energy_fname is not None:
            obs_collections.FormationEnergyFileName = formation_energy_fname
        
        obs_collections.ReadFormationEnergies()

        html_writer.write('<br><b>List of compounds for training</b>')
        html_writer.insert_toggle(start_here=True)
        obs_collections.AddFormationEnergies()
        html_writer.div_end()

        obs_collections.AddRedoxCarriers(anchored=anchored_redox)

        html_writer.write('</br><b>List of NIST reactions for training</b>')
        html_writer.insert_toggle(start_here=True)
        obs_collections.AddNistDatabase()
        html_writer.div_end()
        
        return obs_collections
    
    def AddObservation(self, obs_id, obs_type, url, anchored, dG0, sparse):
        obs = KeggObservation(obs_id, obs_type, url, anchored, dG0, sparse)
        obs.Normalize()
        self.observations.append(obs)

    def ReadFormationEnergies(self):
        """
            Reads the entire table of formation energies which are to be used
            later both to add them directly to the observed data table and to
            be used for normalizing NIST data.
        """

        self.formation_dict = {}
        
        for label in ['training', 'testing']:
            ptable = PsuedoisomerTableThermodynamics.FromCsvFile(self.FormationEnergyFileName,
                                                                 label=label)
            for cid in ptable.get_all_cids():
                pmatrix = ptable.cid2PseudoisomerMap(cid).ToMatrix() 
                if len(pmatrix) != 1:
                    raise Exception("multiple training species for C%05d" % cid)
                nH, charge, nMg, dG0 = pmatrix[0]
                if cid in self.cid2nH_nMg:
                    if (nH, nMg) != self.cid2nH_nMg[cid]:
                        raise Exception("The pseudoisomer of C%05d "
                                        "in the formation energy table (nH=%d) "
                                        "is not consistent with the pKa table (nH=%d)."
                                        % (cid, nH, self.cid2nH_nMg[cid][0]))
                else:
                    self.cid2nH_nMg[cid] = (nH, nMg) 
                diss_table = self.dissociation.GetDissociationTable(cid, False)
                if diss_table is None:
                    raise Exception("C%05d has no pKa data" % cid)
                diss_table.SetFormationEnergyByNumHydrogens(dG0, nH, nMg)
                dG0_prime = diss_table.Transform(pH=self.pH, I=self.I, pMg=self.pMg, T=self.T)
                ref = ptable.cid2SourceString(cid)
                self.formation_dict[cid] = (label, ref, dG0_prime,
                                            dG0, nH, charge, nMg)
    
    def AddFormationEnergies(self):
        """
            Add observations based on a table of derived chemical formation energies.
            If working in self.transformed mode, the dG0_f is Legendre-transformed
            to dG0'_f using the table of dissociation constants.
        """
        
        for cid, v in self.formation_dict.iteritems():
            label, ref, dG0_prime, dG0, nH, _charge, nMg = v
            obs_id = "%s (%d)" % (self.kegg.cid2name(cid), nH)
            logging.debug('Adding the formation energy of %s', obs_id)
            sparse = {cid:1}
    
            html_text = ""
            html_text += "<b>%s [%s]</b></br>\n" % (obs_id, ref)
            html_text += '<font size="1">\n'
            html_text += 'KEGG ID = C%05d</br>\n' % cid
            if self.transformed:
                html_text += "%s = %.1f, " % (self.gibbs_symbol, dG0_prime)
                html_text += 'pH = %g, I = %g, pMg = %g, T = %g</br>\n' % \
                    (self.pH, self.I, self.pMg, self.T)
            else:
                html_text += "%s = %.1f, " % (self.gibbs_symbol, dG0)
                html_text += 'nH = %d, nMg = %d</br>\n' % (nH, nMg)
            html_text += 'Stoichiometry = %s</br>\n' % str(sparse)

            if self.transformed:
                self.AddObservation(obs_id=obs_id, obs_type=KeggObservation.TYPE_FORMATION,
                    anchored=(label=='testing'), dG0=dG0_prime, sparse=sparse,
                    url=self.kegg.cid2link(cid))
            else:
                self.AddObservation(obs_id=obs_id, obs_type=KeggObservation.TYPE_FORMATION,
                    anchored=(label=='testing'), dG0=dG0, sparse=sparse,
                    url=self.kegg.cid2link(cid))
            html_text += '</font>\n'
            self.html_writer.write(html_text)

    def AddNistDatabase(self, assert_decomposition=True):
        """
            Add the observations based on equilibrium constants from the NIST database.
            If using non-transformed group contribution, it is required to reverse Legendre-
            transform the data in order to get chemical reaction energies.
            
            This methods tries to use the same pseudoisomers for each compound as
            the table of formation energies, in order to have less columns in the 
            final stoichiometric matrix.
        """
        nist = Nist()
        
        # create a dictionary from each unique reaction to the list of measured dG0'
        # and subtract from dG0' the formation energies of the anchored compounds
        for r, nist_row_data in enumerate(nist.SelectRowsFromNist()):
            obs_id = "NIST%03d" % r
            msg = ""
            try:
                if not self.transformed:
                    dG0 = self.dissociation.ReverseTransformNistRow(
                                                nist_row_data, self.cid2nH_nMg)
                else:
                    dG0 = nist_row_data.dG0_r # we are using transformed energies

                self.AddObservation(obs_id=obs_id, obs_type=KeggObservation.TYPE_REACTION,
                                    anchored=False, dG0=dG0,
                                    sparse=nist_row_data.reaction.sparse,
                                    url=nist_row_data.url)

            except MissingDissociationConstantError as e:
                msg = 'Cannot reverse transform NIST%03d because of' \
                      ' of a missing dissociation constant for C%05d' % (r, e.cid)
                logging.warning(msg)
                dG0 = None
            
            html_text = ""
            html_text += "<b id=%s>%s</b></br>\n" % (obs_id, obs_id)
            html_text += '<font size="1">\n'
            html_text += "NIST conditions: pH = %g, I = %g, pMg = %g, T = %g</br>\n" % \
                         (nist_row_data.pH, nist_row_data.I,
                          nist_row_data.pMg, nist_row_data.T)
            html_text += 'NIST reference: <a href="%s">%s</a></br>\n' % \
                         (nist_row_data.url, nist_row_data.ref_id)
            html_text += 'EC = %s</br>\n' % nist_row_data.ec
            html_text += "Reaction: %s</br>\n" % \
                         nist_row_data.reaction.to_hypertext(show_cids=False)
            html_text += "%s = %.1f</br>\n" % \
                (thermodynamic_constants.symbol_dr_G0_prime,
                 nist_row_data.dG0_r)
            
            if dG0 is None:
                html_text += 'WARNING: %s</br>\n' % msg
            else:
                html_text += '%s = %.1f</br>\n' % \
                             (self.gibbs_symbol, dG0)
            html_text += '</font></br></br>\n'
            self.html_writer.write(html_text)

    def AddRedoxCarriers(self, anchored=True):
        redox_carriers = RedoxCarriers()
        for name, rc in redox_carriers.iteritems():
            # make sure all redox carriers have a pKa table.
            # for some which don't (usually because their structure is not explicit)
            # we assume that the table is empty (i.e. no pKas exist).
            if self.dissociation.GetDissociationTable(rc.cid_ox) is None:
                self.dissociation.SetOnlyPseudoisomer(rc.cid_ox, rc.nH_ox, rc.z_ox)
            if self.dissociation.GetDissociationTable(rc.cid_red) is None:
                self.dissociation.SetOnlyPseudoisomer(rc.cid_red, rc.nH_red, rc.z_red)

            obs_id = name + " redox"
            self.cid2nH_nMg[rc.cid_ox] = (rc.nH_ox, 0)
            self.cid2nH_nMg[rc.cid_red] = (rc.nH_red, 0)
            sparse = {rc.cid_ox:-1,
                      rc.cid_red:1}

            dG0 = rc.ddG0_prime
            if not self.transformed:
                reaction = Reaction(obs_id, sparse)
                
                ddG0 = self.dissociation.ReverseTransformReaction(
                    reaction, pH=rc.pH, I=0, pMg=default_pMg, T=default_T,
                    cid2nH_nMg=self.cid2nH_nMg)
                dG0 -= ddG0
            
            self.AddObservation(obs_id=obs_id, obs_type=KeggObservation.TYPE_REDOX, url="",
                                anchored=anchored, dG0=dG0, sparse=sparse)
            
    def GetStoichiometry(self):
        """ 
            Returns:
                cids         - a list of pseudoisomer IDs
                obs_ids      - a list of observation IDs
                S            - the observation stoichiometric matrix
                gibbs_values - a row vector of the dG0s
                anchored     - a row vector indicating which obs. are anchored 
        """
        n = len(self.observations) # number of observations
        
        cids = set()
        for obs in self.observations:
            cids.update(obs.sparse.keys())
        cids = sorted(cids)
        
        # create the stoichiometric matrix S (rows=pseudoisomers, cols=observations)
        S = np.matrix(np.zeros((len(cids), n)))
        gibbs_values = np.matrix(np.zeros((1, n)))
        anchored = np.matrix(np.zeros((1, n)))
        obs_ids = []
        for i_obs, obs in enumerate(self.observations):
            obs_ids.append(obs.obs_id)
            for cid, coeff in obs.sparse.iteritems():
                i_cid = cids.index(cid)
                S[i_cid, i_obs] = coeff
            gibbs_values[0, i_obs] = obs.dG0
            if obs.anchored:
                anchored[0, i_obs] = 1
        
        return cids, S, gibbs_values, anchored
                
    def ToDatabase(self, db, table_name):
        # the table 'group_observation' will contain all the observed data
        # that is used for training later
        db.CreateTable(table_name,
                      'id TEXT, type TEXT, url TEXT, anchored BOOL, dG0 REAL, reaction TEXT',
                      drop_if_exists=True)

        for obs in self.observations:
            obs.ToDatabase(db, table_name)
        
        db.Commit()
        
    def ToCSV(self, obs_fname):
        # write all observations
        csv_writer = csv.writer(open(obs_fname, 'w'))
        csv_writer.writerow(['id', 'type', 'url', 'anchored', 'dG0', 'reaction'])
        for obs in self.observations:
            csv_writer.writerow(obs.ToCsvRow())
        
    @staticmethod
    def FromDatabase(db, table_name, transformed=False):
        html_writer = NullHtmlWriter()
        dissociation = None
        obs_collections = KeggObervationCollection(html_writer, 
                                                    dissociation, transformed)
        obs_collections.observations = []
        for row in db.DictReader(table_name):
            obs = KeggObservation.FromDatabaseRow(row)
            obs_collections.observations.append(obs)
        return obs_collections
        
    def ReportToHTML(self):
        rowdicts = []
        for i, obs in enumerate(self.observations):
            rowdict = {'#': i, 'ID': obs.obs_id,
                       'type': obs.obs_type, 'anchored': obs.anchored}
            rowdict[self.gibbs_symbol] = '%.1f' % obs.dG0
            rowdict['Reaction'] = str(self.observations[i].sparse)
            rowdicts.append(rowdict)

        self.html_writer.write('</br><b>Observation summary</b>')
        self.html_writer.insert_toggle(start_here=True)
        self.html_writer.write('<font size="1">\n')
        self.html_writer.write_table(rowdicts, 
                        headers=['#', 'ID', 'type', 'anchored', 
                                 'Reaction', self.gibbs_symbol])
        self.html_writer.write('</font>\n')
        self.html_writer.div_end()

    def ConvertFormation2Reaction(self, output_fname):
        logging.info("Converting all formation energies to reactions")
        output_csv = csv.writer(open(output_fname, 'w'))
        
        # keep the format used for TECRDB
        output_csv.writerow(('ref', 'ID', 'method', 'eval', 'EC', 'name',
                             'kegg_reaction', 'reaction',
                             'dG0\'', 'T', 'I', 'pH', 'pMg'))
        
        atom2cid = {}
        for atom, (name, stoich) in KeggObservation.ATOM2ELEMENT.iteritems():
            cid, _, _ = self.kegg.name2cid(name, 0)
            if cid is None:
                raise Exception("Cannot find the element %s in the KEGG database" % name)
            atom2cid[atom] = (cid, stoich)
            #output_csv.writerow(('element',
            #                     'C%05d' % cid, 'formation', 'A', '',
            #                     'formation of %s' % self.kegg.cid2name(cid),
            #                     "C%05d" % cid,
            #                     name, 0, self.T, self.I, self.pH, self.pMg))

        
        for label in ['training', 'testing']:
            ptable = PsuedoisomerTableThermodynamics.FromCsvFile(self.FormationEnergyFileName,
                                                                 label=label)
            for cid in ptable.get_all_cids():
                pmatrix = ptable.cid2PseudoisomerMap(cid).ToMatrix() 
                if len(pmatrix) != 1:
                    raise Exception("multiple training species for C%05d" % cid)
                nH, charge, nMg, dG0 = pmatrix[0]
                diss_table = dissociation.GetDissociationTable(cid, False)
                if diss_table is None:
                    continue
                diss_table.SetFormationEnergyByNumHydrogens(dG0, nH, nMg)
                dG0_prime = diss_table.Transform(pH=self.pH, I=self.I, pMg=self.pMg, T=self.T)
                ref = ptable.cid2SourceString(cid)
                
                atom_bag = self.kegg.cid2atom_bag(cid)
                if not atom_bag:
                    continue
                
                ne = self.kegg.cid2num_electrons(cid)
                elem_ne = 0
                sparse = {cid: 1}
                for elem, count in atom_bag.iteritems():
                    if elem == 'H':
                        continue
                    elem_ne += count * Molecule.GetAtomicNum(elem)
                    elem_cid, elem_coeff = atom2cid[elem]
                    sparse.setdefault(elem_cid, 0)
                    sparse[elem_cid] += -count * elem_coeff
                
                # use the H element to balance the electrons in the formation
                # reactions (we don't need to balance protons since this is 
                # a biochemical reaction, so H+ are 'free').
                H_cid, H_coeff = atom2cid['H']
                sparse[H_cid] = (elem_ne - ne) * H_coeff
                reaction = Reaction("formation of %s" % self.kegg.cid2name(cid), sparse)
                
                output_csv.writerow((ref,
                                     'C%05d' % cid, 'formation', 'A', '',
                                     'formation of %s' % self.kegg.cid2name(cid),
                                     reaction.FullReactionString(),
                                     reaction.FullReactionString(show_cids=False),
                                     '%.2f' % dG0_prime, self.T, self.I, self.pH, self.pMg))


if __name__ == "__main__":
    dissociation = DissociationConstants.FromPublicDB()
    html_writer = NullHtmlWriter()
    obs_col = KeggObervationCollection(html_writer, dissociation, transformed=True)
    
    obs_col.ConvertFormation2Reaction(output_fname='../res/formation_reactions.csv')