import unittest
from pygibbs.groups_data import GroupsData
from pygibbs import group_decomposition
from toolbox.molecule import Molecule
import logging


PHOSPHATE = Molecule.FromSmiles('[O-]P([O-])(=O)O')
ATP = Molecule.FromSmiles('C1=NC2=C(C(=N1)N)N=CN2C3C(C(C(O3)COP(=O)([O-])OP(=O)([O-])OP(=O)([O-])O)O)O')
A4P = Molecule.FromSmiles('C1=NC2=C(C(=N1)N)N=CN2C3C(C(C(O3)COP(=O)([O-])OP(=O)([O-])OP(=O)([O-])OP(=O)([O-])O)O)O')

class GroupsDecompositionTest(unittest.TestCase):
    """Tests for GroupsDecomposition"""
    
    def setUp(self):
        self.groups_decomposer = group_decomposition.GroupDecomposer.FromGroupsFile(
            open('../data/thermodynamics/groups_species.csv', 'r'))
    
    def testFindPhosphateChains(self):
        ps = group_decomposition.GroupDecomposer.FindPhosphateChains(PHOSPHATE,
                                                                     ignore_protonations=False)
        
        for unused_grp, l in ps:
            self.assertTrue(not l)
        
        mk_ps_dict = lambda ps: dict((key, l) for key, l in ps)
        mk_ps_string = lambda ps: ', '.join(["%s x %d" % (str(key), len(l)) for key, l in ps if l != []])

        # Find chains in A4P        
        ps = group_decomposition.GroupDecomposer.FindPhosphateChains(A4P,
                                                                     ignore_protonations=False,
                                                                     max_length=4)
        logging.info("A4P: " + mk_ps_string(ps))
        ps_dict = mk_ps_dict(ps)
        
        A4P_phosphate_groups = ['middle H0', 'middle chain H0', 'final H1']
        for group_name in A4P_phosphate_groups:
            l = ps_dict.pop(GroupsData.PHOSPHATE_DICT[group_name])
            self.assertEqual(1, len(l))
        for l in ps_dict.itervalues():
            self.assertFalse(l)
        
        # Find chains in ATP
        ps = group_decomposition.GroupDecomposer.FindPhosphateChains(ATP,
                                                                     ignore_protonations=False,
                                                                     max_length=4)
        logging.info("ATP: " + mk_ps_string(ps))
        ps_dict = mk_ps_dict(ps)
        A4P_phosphate_groups = ['middle chain H0', 'final H1']
        for group_name in A4P_phosphate_groups:
            v = ps_dict.pop(GroupsData.PHOSPHATE_DICT[group_name])
            self.assertEqual(1, len(v))
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
    
def Suite():
    suites = (unittest.makeSuite(GroupsDecompositionTest,'test'))
    return unittest.TestSuite(suites)

    
if __name__ == '__main__':
    unittest.main()