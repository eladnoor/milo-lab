import unittest
from pygibbs import groups_data
from pygibbs import group_decomposition
from toolbox.molecule import Molecule


PHOSPHATE = Molecule.FromSmiles('OP(O)(=O)O')
ATP = Molecule.FromSmiles('C1=NC2=C(C(=N1)N)N=CN2C3C(C(C(O3)COP(=O)(O)OP(=O)(O)OP(=O)(O)O)O)O')
A4P = Molecule.FromSmiles('C1=NC2=C(C(=N1)N)N=CN2C3C(C(C(O3)COP(=O)(O)OP(=O)(O)OP(=O)(O)OP(=O)(O)O)O)O')

class GroupsDecompositionTest(unittest.TestCase):
    """Tests for GroupsDecomposition"""
    
    def setUp(self):
        self.groups_decomposer = group_decomposition.GroupDecomposer.FromGroupsFile(
            '../data/thermodynamics/groups_species.csv')
    
    def testFindPhosphateChains(self):
        ps = group_decomposition.GroupDecomposer.FindPhosphateChains(PHOSPHATE,
                                                                     ignore_protonations=True)
        
        for unused_grp, l in ps:
            self.assertTrue(not l)
        
        mk_ps_dict = lambda ps: dict((key, l) for key, l in ps)

        # Find chains in A4P        
        ps = group_decomposition.GroupDecomposer.FindPhosphateChains(A4P,
                                                                     ignore_protonations=True,
                                                                     max_length=4)
        ps_dict = mk_ps_dict(ps)
        
        self.assertEqual(1, len(ps_dict[groups_data.GroupsData.MIDDLE_P_2]))
        self.assertEqual(1, len(ps_dict[groups_data.GroupsData.MIDDLE_2_PHOSPHATE]))
        self.assertEqual(1, len(ps_dict[groups_data.GroupsData.FINAL_P_1]))
        ps_dict.pop(groups_data.GroupsData.MIDDLE_P_2)
        ps_dict.pop(groups_data.GroupsData.MIDDLE_2_PHOSPHATE)
        ps_dict.pop(groups_data.GroupsData.FINAL_P_1)
        for l in ps_dict.itervalues():
            self.assertFalse(l)
        
        # Find chains in ATP
        ps = group_decomposition.GroupDecomposer.FindPhosphateChains(ATP,
                                                                     ignore_protonations=True,
                                                                     max_length=4)
        ps_dict = mk_ps_dict(ps)
        self.assertEqual(1, len(ps_dict[groups_data.GroupsData.MIDDLE_2_PHOSPHATE]))
        self.assertEqual(1, len(ps_dict[groups_data.GroupsData.FINAL_P_1]))
        ps_dict.pop(groups_data.GroupsData.MIDDLE_2_PHOSPHATE)
        ps_dict.pop(groups_data.GroupsData.FINAL_P_1)
        for l in ps_dict.itervalues():
            self.assertFalse(l)    
    
    def testDecompose(self):
        # Phosphate has no decomposition.
        decomposition = self.groups_decomposer.Decompose(PHOSPHATE,
                                                         ignore_protonations=True)

        for grp in decomposition.groups:
            self.assertTrue(not grp[-1])
        
        # When we're strict, should lead to an exception.
        self.assertRaises(group_decomposition.GroupDecompositionError,
                          self.groups_decomposer.Decompose, PHOSPHATE,
                          ignore_protonations=True, strict=True)
        
    def testPseudoisomerVectors(self):
        # Phosphate has no decomposition, hence no pseudoisomers.
        decomposition = self.groups_decomposer.Decompose(PHOSPHATE,
                                                         ignore_protonations=True)
        pseudoisomers = decomposition.PseudoisomerVectors()
        self.assertFalse(pseudoisomers)
    
    
if __name__ == '__main__':
    unittest.main()