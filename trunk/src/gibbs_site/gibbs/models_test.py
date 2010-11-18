#!/usr/bin/python

import unittest
from util import django_utils

# NOTE(flamholz): This is crappy. We're using the real database for
# a unit test. I wish I knew of a better way.
django_utils.SetupDjango()

from gibbs import models

class SpeciesFormationEnergyTest(unittest.TestCase):
    """Tests for SpeciesFormationEnergy"""
        
    def testTransform(self):
        specie = models.Specie(number_of_hydrogens=2, net_charge=-1,
                               formation_energy=-10.0)
        
        # Test some hand-calculated numbers: (ph, ionic strength, result)
        test_data = (
             # Move pH, keep ionic strength constant.
             (6.0, 0.1, 59.071),
             (6.5, 0.1, 64.776),
             (7.0, 0.1, 70.481),
             (7.5, 0.1, 76.186),
             (8.0, 0.1, 81.891),
             # Move ionic strength, keep pH constant.
             (7.0, 0.0001, 69.898),
             (7.0, 0.001, 69.957),
             (7.0, 0.01, 70.121),
             (7.0, 0.05, 70.349),
             (7.0, 0.11, 70.501),
             (7.0, 0.15, 70.566),
             (7.0, 0.2, 70.629),
             # Move both.
             (6.0, 0.0001, 58.488),
             (6.5, 0.001, 64.252),
             (8.0, 0.0001, 81.308),
             (7.5, 0.15, 76.271),
             )

        for ph, ionic_strength, expected_transform in test_data:
            actual_transform = specie.Transform(pH=ph, ionic_strength=ionic_strength)
            self.assertAlmostEqual(expected_transform, actual_transform, 3)


class CompoundTest(unittest.TestCase):
    
    def testGetAtomBag(self):
        compound = models.Compound(kegg_id='fake compound')
        self.assertEqual(None, compound.GetAtomBag())
        
        compound.formula = 'C12H22O11'
        expected_atom_bag = {'C': 12, 'H': 22, 'O': 11}
        self.assertEqual(expected_atom_bag, compound.GetAtomBag())
        
        compound.formula = 'C10H16N5O12P3S'
        expected_atom_bag = {'C': 10, 'H': 16, 'N': 5,
                             'O': 12, 'P': 3, 'S': 1}
        self.assertEqual(expected_atom_bag, compound.GetAtomBag())

        # Contains an R group.
        compound.formula = 'C10R11H16N5O12P3S'
        self.assertEqual(None, compound.GetAtomBag())

    def testDeltaG(self):
        # Create a test compound.
        species = [models.Specie(number_of_hydrogens=12, net_charge=0,
                                 formation_energy=-10.5),
                   models.Specie(number_of_hydrogens=11, net_charge=-1,
                                 formation_energy=-12.1),
                   models.Specie(number_of_hydrogens=10, net_charge=-2,
                                 formation_energy=-13.4)]

        compound = models.Compound(kegg_id='fake compound')
        compound._all_species = species
        
        # Format: ph, ionic strength, dG.
        test_data = ((6.5, 0.1, 361.094),
                     (7.0, 0.1, 389.619),
                     (7.5, 0.1, 418.143),
                     (7.0, 0.0001, 386.118),
                     (7.0, 0.001, 386.473),
                     (7.0, 0.2, 390.505))

        for ph, i_s, expected_dg in test_data:
            self.assertAlmostEqual(expected_dg,
                                   compound.DeltaG(pH=ph, ionic_strength=i_s),
                                   3)


if __name__ == '__main__':
    unittest.main()