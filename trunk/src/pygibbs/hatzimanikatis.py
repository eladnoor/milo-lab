import csv, sys
import numpy as np
from matplotlib import pyplot

from pygibbs.thermodynamics import Thermodynamics, MissingCompoundFormationEnergy
from pygibbs.thermodynamic_constants import J_per_cal
from pygibbs.dissociation_constants import DissociationConstants
from pygibbs.kegg import Kegg
from pygibbs.thermodynamic_constants import R
from toolbox.database import SqliteDatabase
from pygibbs.pseudoisomer import PseudoisomerMap
from toolbox.util import _mkdir
from pygibbs.kegg_reaction import Reaction
from toolbox.linear_regression import LinearRegression

HATZI_CSV_FNAME = "../data/thermodynamics/hatzimanikatis_cid.csv"

class Hatzi (Thermodynamics):
    
    def __init__(self, use_pKa=True):
        if use_pKa:
            Thermodynamics.__init__(self, "Jankowski et al. (+pKa)")
            self.dissociation = DissociationConstants.FromPublicDB()
        else:
            Thermodynamics.__init__(self, "Jankowski et al.")
            self.dissociation = None
        self.db = SqliteDatabase('../res/gibbs.sqlite', 'w')
        self.cid2pmap_dict = {}
        
        # the conditions in which Hatzimanikatis makes his predictions
        self.Hatzi_pH = 7.0
        self.Hatzi_I = 0.0
        self.Hatzi_pMg = 14.0
        self.Hatzi_T = 298.15
        
        self.kegg = Kegg.getInstance()

        # for some reason, Hatzimanikatis doesn't indicate that H+ is zero,
        # so we add it here
        H_pmap = PseudoisomerMap()
        H_pmap.Add(0, 0, 0, 0)
        self.SetPseudoisomerMap(80, H_pmap)

        self.cid2dG0_tag_dict = {80 : 0}
        self.cid2charge_dict = {80 : 0}

        for row in csv.DictReader(open(HATZI_CSV_FNAME, 'r')):
            cid = int(row['ENTRY'][1:])
            self.cid2source_string[cid] = 'Jankowski et al. 2008'
            if row['DELTAG'] == "Not calculated":
                continue
            if cid == 3178:
                # this compound, which is supposed to be "Tetrahydroxypteridine"
                # seems to be mapped to something else by Hatzimanikatis
                continue
            self.cid2dG0_tag_dict[cid] = float(row['DELTAG']) * J_per_cal
            self.cid2charge_dict[cid] = int(row['CHARGE'])

    def charge2nH(self, cid, charge):
        nH_z_pair = self.kegg.cid2nH_and_charge(cid)
        if nH_z_pair:
            nH, z = nH_z_pair
            return nH + (charge - z)
        else:
            raise ValueError("cannot infer the nH of C%05d" % cid)

    def GeneratePseudoisomerMap(self, cid):
        """
            Generates a PseudoisomerMap for the CID using the
            data in Hatzimanikatis' CID table. Note that when
            using the dissociation constant table to generate all
            psuedoisomers, one must reverse-transform the dG0
            since Hatzimanikatis is using dG0' and not dG0. It is 
            misleading because in his model H+ does have a dG0'
            associated to it, according to the pH=7, and that must be
            subtracted 
        """ 
        
        charge = self.cid2charge_dict[cid]
        dG0_tag = self.cid2dG0_tag_dict[cid]
        if self.dissociation is not None:
            diss_table = self.dissociation.GetDissociationTable(cid)
            if diss_table.min_nH == None:
                try:
                    nH = self.charge2nH(cid, charge)
                except KeyError:
                    raise MissingCompoundFormationEnergy('The compound C%05d is missing from KEGG' % cid)
                except ValueError:
                    nH = 0
                diss_table.min_nH = nH 
                diss_table.min_charge = charge
            else:
                nH = diss_table.min_nH + (charge - diss_table.min_charge)
            
            dG0_hplus = -R * self.T * np.log(10) * self.Hatzi_pH
            dG0_tag -= nH*dG0_hplus
            
            diss_table.SetTransformedFormationEnergy(
                dG0_tag, pH=self.Hatzi_pH, I=self.Hatzi_I, 
                pMg=self.Hatzi_pMg, T=self.Hatzi_T)
            
            return diss_table.GetPseudoisomerMap()
        else:
            pmap = PseudoisomerMap()
            pmap.Add(nH=charge, z=charge, nMg=0, dG0=dG0_tag)
            return pmap

    def SetPseudoisomerMap(self, cid, pmap):
        self.cid2pmap_dict[cid] = pmap

    def AddPseudoisomer(self, cid, nH, z, nMg, dG0):
        self.cid2pmap_dict.setdefault(cid, PseudoisomerMap())
        self.cid2pmap_dict[cid].Add(nH, z, nMg, dG0)

    def cid2PseudoisomerMap(self, cid):
        if cid in self.cid2pmap_dict:
            return self.cid2pmap_dict[cid]
        if cid in self.cid2dG0_tag_dict:
            pmap = self.GeneratePseudoisomerMap(cid)
            self.SetPseudoisomerMap(cid, pmap)
            return pmap

        raise MissingCompoundFormationEnergy(
            "The compound C%05d does not have a value for its"
            " formation energy of any of its pseudoisomers" % cid, cid)

    def get_all_cids(self):
        return sorted(self.cid2dG0_tag_dict.keys())
        

def TestGroupMatrix():
    group_filename = '../data/thermodynamics/hatzimanikatis_groups.csv'
    all_group_names = []
    sparse_matrix = []
    dG_vector = []
    line_no = 0
    for row in csv.DictReader(open(group_filename)):
        line_no += 1
        if row['est_dG'] == "None":
            continue
        dG_vector.append(float(row['est_dG']))

        sparse_groupvec = []
        if row['groups'] != "":
            for token in row['groups'].split(' | '):
                try:
                    [group_name, coeff] = token.split(' : ', 1)
                except ValueError:
                    raise Exception("cannot parse this token (line %d): %s\n" % 
                                    (line_no, token))
                coeff = float(coeff)
                if group_name not in all_group_names:
                    all_group_names.append(group_name)
                group_index = all_group_names.index(group_name)
                sparse_groupvec.append((group_index, coeff))
        sparse_matrix.append(sparse_groupvec)
        
    full_matrix = np.zeros((len(sparse_matrix), len(all_group_names)))
    for i in range(len(sparse_matrix)):
        for j, coeff in sparse_matrix[i]:
            full_matrix[i, j] = coeff
    dG_vector = np.array(dG_vector, ndmin=2).T
    
    #print full_matrix.shape
    #print LinearRegression.MatrixRank(full_matrix)
    #print dG_vector.shape
    
    #augmented_matrix = np.hstack([full_matrix, dG_vector])
    #_U, s, _V = np.linalg.svd(augmented_matrix, full_matrices=False)
    #print sorted(s)
    
    contributions, _K = LinearRegression.LeastSquares(full_matrix, dG_vector)
    for i, group_name in enumerate(all_group_names):
        print "%s,%.3f" % (group_name, contributions[i, 0])
        
    pyplot.plot(dG_vector, dG_vector-np.dot(full_matrix, contributions), '.')
    pyplot.show()

if __name__ == "__main__":
    if len(sys.argv) > 1 and sys.argv[1] == 'test':
        TestGroupMatrix()
        sys.exit(0)
    
    _mkdir('../res')
    db = SqliteDatabase('../res/gibbs.sqlite', 'w')
    H_nopka = Hatzi(use_pKa=False)
    H_withpka = Hatzi(use_pKa=True)
    H_withpka.ToDatabase(db, 'hatzi_thermodynamics')
    
    #H.ToDatabase(db, 'hatzi_gc')
    #H.I = 0.25
    #H.T = 300;
    #sparse_reaction = {13:-1, 1:-1, 9:2}
    #sparse_reaction = {36:-1, 3981:1}
    #sparse_reaction = {6:-1, 143:-1, 234:1, 5:1}
    #sparse_reaction = {1:-1, 499:-1, 603:1, 86:1}
    #sparse_reaction = {1:-1, 6:-1, 311:-1, 288:1, 5:1, 80:2, 26:1}
    #sparse_reaction = {408:-1, 6:-1, 4092:1, 5:1}
    #sparse_reaction = {588:-1, 1:-1, 114:1, 9:1}
    #sparse_reaction = {1:-1, 3:-1, 149:-1, 288:1, 4:1, 80:2, 22:1}
    react = Reaction("reaction", {408:-1, 6:-1, 4092:1, 5:1})
    
    #sys.stdout.write("The dG0_r of PPi + H20 <=> 2 Pi: \n\n")
    
    react.Balance()
    print react.FullReactionString()
    
    sys.stdout.write("%5s | %5s | %6s | %6s\n" % ("pH", "I", "T", "dG0_r"))
    for pH in np.arange(5, 10.01, 0.25):
        H_withpka.pH = pH
        sys.stdout.write("%5.2f | %5.2f | %6.1f | %6.2f\n" % 
                         (H_withpka.pH, H_withpka.I, H_withpka.T, 
                          react.PredictReactionEnergy(H_withpka)))

    for cid in react.get_cids():
        print '-'*50
        print "C%05d - %s:" % (cid, H_withpka.kegg.cid2name(cid))
        print H_withpka.kegg.cid2inchi(cid)
        print "Pseudoisomers:\n", H_withpka.cid2PseudoisomerMap(cid)
        print "dG0'_f = %.1f kJ/mol" % H_withpka.cid2PseudoisomerMap(cid).Transform(pH=7, I=0, pMg=10, T=298.15)
        print "dG0_f = %.1f kJ/mol" % H_nopka.cid2PseudoisomerMap(cid).Transform(pH=7, I=0, pMg=10, T=298.15)
        print "Dissociations:\n", H_withpka.dissociation.GetDissociationTable(cid)
