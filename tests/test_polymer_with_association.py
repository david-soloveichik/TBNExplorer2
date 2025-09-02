"""Tests for polymer free energy including association penalty."""

import numpy as np
import pytest

from tbnexplorer2.model import TBN, BindingSite, Monomer
from tbnexplorer2.polymer_basis import Polymer, compute_assoc_energy_penalty


class TestPolymerWithAssociation:
    def test_free_energy_with_association_penalty(self):
        """Test free energy computation with association penalty."""
        # Create TBN with two monomers: a and a*
        monomer1 = Monomer(name="M1", binding_sites=[BindingSite("a", False)], concentration=100, original_line="M1: a")
        monomer2 = Monomer(name="M2", binding_sites=[BindingSite("a", True)], concentration=100, original_line="M2: a*")
        tbn = TBN([monomer1, monomer2], {"a": 0})

        # Create dimer polymer (one of each monomer)
        polymer = Polymer(np.array([1, 1]), [monomer1, monomer2], tbn)

        # Test with non-zero association parameters
        deltaG = [5.0, 3.0]  # dG_assoc=5, dH_assoc=3
        temperature = 37.0

        # Association penalty for 2 monomers = bimolecular * (2-1) = bimolecular
        assoc_penalty = compute_assoc_energy_penalty(2, temperature, 5.0, 3.0)

        expected_free_energy = assoc_penalty
        actual_free_energy = polymer.compute_free_energy(deltaG, temperature)

        assert abs(actual_free_energy - expected_free_energy) < 1e-6

    def test_free_energy_larger_polymer_with_association(self):
        """Test free energy for larger polymer with association penalty."""
        # Create TBN with monomers that can form larger polymers
        monomer1 = Monomer(
            name="M1",
            binding_sites=[BindingSite("a", False), BindingSite("b", False)],
            concentration=100,
            original_line="M1: a b",
        )
        monomer2 = Monomer(
            name="M2",
            binding_sites=[BindingSite("a", True), BindingSite("b", True)],
            concentration=100,
            original_line="M2: a* b*",
        )
        tbn = TBN([monomer1, monomer2], {"a": 0, "b": 1})

        # Create polymer with 3 M1 and 3 M2 (total 6 monomers)
        polymer = Polymer(np.array([3, 3]), [monomer1, monomer2], tbn)

        deltaG = [4.0, 2.0]  # dG_assoc=4, dH_assoc=2
        temperature = 25.0

        # Association penalty for 6 monomers
        assoc_penalty = compute_assoc_energy_penalty(6, temperature, 4.0, 2.0)

        expected_free_energy = assoc_penalty
        actual_free_energy = polymer.compute_free_energy(deltaG, temperature)

        assert abs(actual_free_energy - expected_free_energy) < 1e-6

    def test_free_energy_temperature_dependence(self):
        """Test that free energy changes with temperature."""
        # Create simple dimer
        monomer1 = Monomer(name="M1", binding_sites=[BindingSite("a", False)], concentration=100, original_line="M1: a")
        monomer2 = Monomer(name="M2", binding_sites=[BindingSite("a", True)], concentration=100, original_line="M2: a*")
        tbn = TBN([monomer1, monomer2], {"a": 0})
        polymer = Polymer(np.array([1, 1]), [monomer1, monomer2], tbn)

        deltaG = [5.0, 3.0]

        # Compute at different temperatures
        energy_25 = polymer.compute_free_energy(deltaG, 25.0)
        energy_37 = polymer.compute_free_energy(deltaG, 37.0)
        energy_50 = polymer.compute_free_energy(deltaG, 50.0)

        # Free energies should be different at different temperatures
        assert energy_25 != energy_37
        assert energy_37 != energy_50

    def test_singleton_no_association_penalty(self):
        """Test that singleton has no association penalty."""
        # Create a singleton polymer
        monomer = Monomer(
            name="M1",
            binding_sites=[BindingSite("a", False), BindingSite("a", True)],
            concentration=100,
            original_line="M1: a a*",
        )
        tbn = TBN([monomer], {"a": 0})
        polymer = Polymer(np.array([1]), [monomer], tbn)

        deltaG = [5.0, 3.0]  # Non-zero association parameters
        temperature = 37.0

        # Singleton has no bonds and no association penalty
        # (total_monomers = 1, so penalty = bimolecular * (1-1) = 0)
        assert polymer.compute_free_energy(deltaG, temperature) == 0.0

    def test_default_no_association_penalty(self):
        """Test that default (no params) results in no association penalty."""
        # Create simple dimer
        monomer1 = Monomer(name="M1", binding_sites=[BindingSite("a", False)], concentration=100, original_line="M1: a")
        monomer2 = Monomer(name="M2", binding_sites=[BindingSite("a", True)], concentration=100, original_line="M2: a*")
        tbn = TBN([monomer1, monomer2], {"a": 0})
        polymer = Polymer(np.array([1, 1]), [monomer1, monomer2], tbn)

        # With default (None), bond term ignored and no association penalty
        assert polymer.compute_free_energy() == 0.0

        # Test with explicit association params = [0.0, 0.0]; includes water density term
        energy_with_assoc = polymer.compute_free_energy([0.0, 0.0])

        # This should be different from 0.0 due to water density contribution and negative
        assert energy_with_assoc != 0.0
        assert energy_with_assoc < 0.0


if __name__ == "__main__":
    pytest.main([__file__])
