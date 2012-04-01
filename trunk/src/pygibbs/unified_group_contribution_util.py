from pygibbs.groups_data import GroupsData
from pygibbs.group_decomposition import GroupDecomposer
from toolbox.molecule import Molecule
import numpy as np
from pygibbs.thermodynamic_constants import default_T, R

class InChI2FormationEnergy(object):
    
    def __init__(self):
        self.transformed = False
        self.RT = R * default_T
        pass
    
    def Initialize(self, db):
        from pygibbs.unified_group_contribution import UnifiedGroupContribution

        ugc = UnifiedGroupContribution(db)
        ugc.LoadGroups(FromDatabase=True)
        ugc.LoadObservations(FromDatabase=True)
        ugc.LoadGroupVectors(FromDatabase=True)
        ugc.LoadData(FromDatabase=True)
        ugc.init()
        
        self.groups_data = ugc.groups_data
        self.group_decomposer = ugc.group_decomposer

        result_dict = ugc._GetContributionData(ugc.S.copy(), ugc.cids,
                                               ugc.b.copy(), ugc.anchored)
        
        self.g_pgc = result_dict['group_contributions']
        self.P_L_pgc = result_dict['pgc_conservations']
    
    def ToDatabase(self, db):
        self.groups_data.ToDatabase(db)
        db.SaveNumpyMatrix('ugc_group_contributions', self.g_pgc.T)
        db.SaveSparseNumpyMatrix('ugc_group_nullspace', self.P_L_pgc)
        db.Commit()
    
    def FromDatabase(self, db):
        self.groups_data = GroupsData.FromDatabase(db,
                                                   transformed= self.transformed)
        self.group_decomposer = GroupDecomposer(self.groups_data)
        
        self.g_pgc = db.LoadNumpyMatrix('ugc_group_contributions').T
        self.P_L_pgc = db.LoadSparseNumpyMatrix('ugc_group_nullspace')
        
    def EstimateGroupVector(self, groupvec):
        gv = np.matrix(groupvec.Flatten())
        dG0 = float(self.g_pgc * gv.T)
        ker = self.P_L_pgc * gv.T
        return dG0, ker
    
    def EstimateInChI(self, inchi):
        mol = Molecule.FromInChI(inchi)
        mol.RemoveHydrogens()
        decomposition = self.group_decomposer.Decompose(mol, 
                            ignore_protonations=False, strict=True)

        nH = decomposition.Hydrogens()
        charge = decomposition.NetCharge()
        nMg = decomposition.Magnesiums()
        groupvec = decomposition.AsVector()
        dG0, ker = self.EstimateGroupVector(groupvec)
        return dG0, nH, charge, nMg, ker
    
    def GenerateAllPseudoisomers(self, dG0, nH, charge, nMg, pKas):
        """
            Given the values of the most abundant species at pH 7 and a list
            of pKa values, generates the data of all other pseudoisomers
        """
        pKa_higher = [x for x in pKas if 7 < x]
        pKa_lower = [x for x in pKas if 7 > x]
        
        pseudoisomer_list = [{'dG0': round(dG0, 1), 'nH': nH,
                              'charge': charge, 'nMg': nMg}]
        
        ddG0 = 0
        for i, pKa in enumerate(sorted(pKa_higher)):
            ddG0 += self.RT * np.log(10) * pKa
            pseudoisomer_list.append({'dG0': round(dG0 + ddG0, 1), 'nH': nH-1-i,
                                      'charge': charge-1-i, 'nMg': nMg})
        
        ddG0 = 0
        for i, pKa in enumerate(sorted(pKa_lower, reverse=True)):
            ddG0 -= self.RT * np.log(10) * pKa
            pseudoisomer_list.append({'dG0': round(dG0 + ddG0, 1), 'nH': nH+1+i,
                                      'charge': charge+1+i, 'nMg': nMg})
        
        return sorted(pseudoisomer_list, key=lambda x: x['nH'])

if __name__ == "__main__":
    
    from toolbox.database import SqliteDatabase
    inchi2dg = InChI2FormationEnergy()
    db_out = SqliteDatabase('../res/ugc.sqlite', 'w')
    
    if False:
        db = SqliteDatabase('../res/gibbs.sqlite', 'w')
        inchi2dg.Initialize(db)
        inchi2dg.ToDatabase(db_out)
    else:
        inchi2dg.FromDatabase(db_out)
    
    inchi_list = ['InChI=1S/C2H5NO2/c3-1-2(4)5/h1,3H2,(H,4,5)',
                  'InChI=1S/C10H16N5O13P3/c11-8-5-9(13-2-12-8)15(3-14-5)10-7(17)6(16)4(26-10)1-25-30(21,22)28-31(23,24)27-29(18,19)20/h2-4,6-7,10,16-17H,1H2,(H,21,22)(H,23,24)(H2,11,12,13)(H2,18,19,20)/p-3']
    
    for inchi in inchi_list:
        print inchi2dg.EstimateInChI(inchi)
    
    