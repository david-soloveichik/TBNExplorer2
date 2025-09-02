import numpy as np
import pytest

from tbnexplorer2.model import TBN, BindingSite, Monomer
from tbnexplorer2.polymer_basis import Polymer


class TestPolymer:
    def test_free_energy_singleton(self):
        """Test free energy computation for singleton polymer."""
        # Create a simple TBN with one monomer: a a*
        binding_sites = [BindingSite("a", False), BindingSite("a", True)]
        monomer = Monomer(name="M1", binding_sites=binding_sites, concentration=100, original_line="M1: a a*")
        tbn = TBN([monomer], {"a": 0})

        # Create singleton polymer (just the monomer itself)
        polymer = Polymer(np.array([1]), [monomer], tbn)

        # For a singleton a a*, there's no bond formed (self-binding excluded)
        # Free energy should be 0 (no bonds, no association penalty for single monomer)
        # When deltaG is None (default), no association penalty is applied
        assert polymer.compute_free_energy() == 0

    def test_free_energy_dimer(self):
        """Test free energy computation for a dimer."""
        # Create TBN with two monomers: a and a*
        monomer1 = Monomer(name="M1", binding_sites=[BindingSite("a", False)], concentration=100, original_line="M1: a")
        monomer2 = Monomer(name="M2", binding_sites=[BindingSite("a", True)], concentration=100, original_line="M2: a*")
        tbn = TBN([monomer1, monomer2], {"a": 0})

        # Create dimer polymer (one of each monomer)
        polymer = Polymer(np.array([1, 1]), [monomer1, monomer2], tbn)

        # With bond term ignored and no association penalty by default
        assert polymer.compute_free_energy() == 0

    def test_free_energy_complex_polymer(self):
        """Test free energy computation for a more complex polymer."""
        # Create TBN with monomers:
        # M1: a b
        # M2: a* b*
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

        # Create polymer with 2 M1 and 2 M2
        polymer = Polymer(np.array([2, 2]), [monomer1, monomer2], tbn)

        # Matrix A = [[1, -1],   # a row: M1 has a, M2 has a*
        #             [1, -1]]   # b row: M1 has b, M2 has b*
        # |A| = [[1, 1], [1, 1]]
        # x = [2, 2]
        # |A| * x = [1*2 + 1*2, 1*2 + 1*2] = [4, 4], sum = 8
        # A * x = [1*2 + (-1)*2, 1*2 + (-1)*2] = [0, 0], sum = 0
        # Bonds are ignored; default energy is 0 without association penalty
        assert polymer.compute_free_energy() == 0

    def test_free_energy_partial_binding(self):
        """Test free energy for polymer with partial binding."""
        # Create TBN with monomers: {a, b} and {a*, c}
        monomer1 = Monomer(
            name="M1",
            binding_sites=[BindingSite("a", False), BindingSite("b", False)],
            concentration=100,
            original_line="M1: a b",
        )
        monomer2 = Monomer(
            name="M2",
            binding_sites=[BindingSite("a", True), BindingSite("c", False)],
            concentration=100,
            original_line="M2: a* c",
        )
        tbn = TBN([monomer1, monomer2], {"a": 0, "b": 1, "c": 2})

        # Create polymer with 1 M1 and 1 M2
        polymer = Polymer(np.array([1, 1]), [monomer1, monomer2], tbn)

        # Matrix A = [[1, -1],  # a row
        #             [1, 0],   # b row
        #             [0, 1]]   # c row
        # x = [1, 1]
        # |A| * x = [2, 1, 1], sum = 4 (total binding sites)
        # A * x = [0, 1, 1], sum = 2 (excess unstar)
        # Bonds are ignored; default energy is 0 without association penalty
        assert polymer.compute_free_energy() == 0

    def test_free_energy_no_tbn_error(self):
        """Test that computing free energy without TBN raises error."""
        # Create polymer without TBN reference
        monomer = Monomer(name="M1", binding_sites=[BindingSite("a", False)], concentration=100, original_line="M1: a")
        polymer = Polymer(np.array([1]), [monomer], tbn=None)

        with pytest.raises(ValueError, match="Cannot compute free energy without TBN"):
            polymer.compute_free_energy()

    def test_free_energy_repeatable(self):
        """Test that repeated free energy computations return the same value for same params."""
        # Create a simple TBN
        monomer1 = Monomer(name="M1", binding_sites=[BindingSite("a", False)], concentration=100, original_line="M1: a")
        monomer2 = Monomer(name="M2", binding_sites=[BindingSite("a", True)], concentration=100, original_line="M2: a*")
        tbn = TBN([monomer1, monomer2], {"a": 0})

        polymer = Polymer(np.array([1, 1]), [monomer1, monomer2], tbn)

        # Compute twice with same parameters
        e1 = polymer.compute_free_energy()
        e2 = polymer.compute_free_energy()
        assert e1 == e2
