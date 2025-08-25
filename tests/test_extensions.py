"""Tests for the extensions module."""

import tempfile
from pathlib import Path

import numpy as np
import pytest

from extensions.canonical_reactions import CanonicalReactionsComputer, Reaction
from extensions.ibot import IBOTAlgorithm
from tbnexplorer2.model import TBN, BindingSite, Monomer
from tbnexplorer2.parser import TBNParser
from tbnexplorer2.polymer_basis import PolymerBasisComputer


def create_test_tbn():
    """Create a simple TBN for testing."""
    binding_sites = {"a": 0, "b": 1, "c": 2}

    monomers = [
        Monomer(
            binding_sites=[BindingSite("a", False), BindingSite("b", True)],
            concentration=None,
            name="M1",
            original_line="M1: a b*",
        ),
        Monomer(
            binding_sites=[BindingSite("b", False), BindingSite("c", True)],
            concentration=None,
            name="M2",
            original_line="M2: b c*",
        ),
        Monomer(
            binding_sites=[BindingSite("c", False), BindingSite("a", True)],
            concentration=None,
            name="M3",
            original_line="M3: c a*",
        ),
    ]

    return TBN(monomers, binding_sites, concentration_units=None)


class TestCanonicalReactions:
    """Test canonical reactions computation."""

    def test_reaction_balanced(self):
        """Test reaction balance checking."""
        # Balanced reaction: 2A -> B + C (2 reactants, 2 products)
        reaction_vec = np.array([-2, 1, 1, 0])
        reaction = Reaction(reaction_vec)
        assert reaction.is_balanced()

        # Unbalanced reaction: A -> B + C (1 reactant, 2 products)
        reaction_vec = np.array([-1, 1, 1, 0])
        reaction = Reaction(reaction_vec)
        assert not reaction.is_balanced()

    def test_reaction_string_representation(self):
        """Test reaction string formatting."""
        reaction_vec = np.array([-1, -1, 2, 0])
        reaction = Reaction(reaction_vec, polymer_names=["A", "B", "C", "D"])
        assert str(reaction) == "A + B -> 2 C"

    def test_setup_matrices(self):
        """Test B and S matrix setup."""
        tbn = create_test_tbn()
        computer = CanonicalReactionsComputer(tbn)

        # Create mock polymer basis
        polymers = [
            np.array([1, 0, 0]),  # Just M1
            np.array([0, 1, 0]),  # Just M2
            np.array([0, 0, 1]),  # Just M3
            np.array([1, 1, 0]),  # M1 + M2
        ]

        on_target_indices = {0, 1}  # First two polymers are on-target

        computer.setup_matrices(polymers, on_target_indices)

        # Check B matrix dimensions
        assert computer.B_matrix.shape == (3, 4)  # 3 monomers, 4 polymers

        # Check S matrix dimensions
        assert computer.S_matrix.shape == (2, 4)  # 2 off-target, 4 polymers

        # Check S matrix selects off-target polymers
        assert computer.S_matrix[0, 2] == 1  # Selects polymer 2
        assert computer.S_matrix[1, 3] == 1  # Selects polymer 3


class TestIBOTAlgorithm:
    """Test IBOT algorithm."""

    def test_ibot_initialization(self):
        """Test IBOT algorithm initialization."""
        tbn = create_test_tbn()
        polymers = [np.array([1, 0, 0]), np.array([0, 1, 0]), np.array([0, 0, 1])]
        on_target_indices = {0}
        reactions = []

        ibot = IBOTAlgorithm(tbn, polymers, on_target_indices, reactions)

        # Check initial concentration exponents
        assert ibot.mu[0] == 1.0  # On-target has μ = 1
        assert ibot.mu[1] == 0.0  # Off-target starts at 0
        assert ibot.mu[2] == 0.0  # Off-target starts at 0

        # Check unassigned tracking
        assert ibot.unassigned_off_target == {1, 2}

    def test_reaction_metrics(self):
        """Test reaction metrics computation."""
        tbn = create_test_tbn()
        polymers = [np.array([1, 0, 0]), np.array([0, 1, 0]), np.array([0, 0, 1])]
        on_target_indices = {0}
        reactions = []

        ibot = IBOTAlgorithm(tbn, polymers, on_target_indices, reactions)

        # Create a reaction: P0 -> P1 + P2
        reaction_vec = np.array([-1, 1, 1])
        reaction = Reaction(reaction_vec)

        metrics = ibot.compute_reaction_metrics(reaction)

        # Both P1 and P2 are unassigned off-target
        assert metrics.novelty == 2
        # Imbalance = 1*μ(P0) - 1*μ(P1) - 1*μ(P2) = 1*1 - 0 - 0 = 1
        assert metrics.imbalance == 1.0
        # Ratio = 1/2 = 0.5
        assert metrics.ratio == 0.5

    def test_ibot_tbn_generation_with_units(self):
        """Test .tbn generation with proper unit conversion."""
        tbn = create_test_tbn()
        polymers = [np.array([1, 0, 0]), np.array([0, 1, 0])]
        on_target_indices = {0, 1}
        reactions = []  # No reactions for simplicity

        ibot = IBOTAlgorithm(tbn, polymers, on_target_indices, reactions)

        # All on-target polymers should have μ = 1
        concentration_exponents = ibot.run()

        # Test .tbn generation with unit conversion
        with tempfile.NamedTemporaryFile(suffix=".tbn", delete=False) as f:
            output_path = Path(f.name)

        try:
            # Generate with c=100 nM
            ibot.generate_tbn_output(output_path, 100, "nM")

            # Read and verify the file
            content = output_path.read_text()
            assert "\\UNITS: nM" in content

            # Verify concentrations are calculated correctly
            # c' = 100 nM = 100e-9 M = 1e-7 M
            # Each monomer should have concentration = count * (1e-7)^1
            # Converted back to nM = 100 nM for each monomer appearance
            lines = content.strip().split("\n")
            monomer_count = 0
            for line in lines:
                if ", " in line and not line.startswith("#"):
                    monomer_count += 1
                    # Extract concentration value
                    conc_str = line.split(", ")[1]
                    conc_val = float(conc_str)
                    # Monomers M1 and M2 each appear once in their respective polymers with μ=1
                    # So concentration should be 100 nM
                    # M3 doesn't appear in any of our polymers, so should be 0
                    if monomer_count <= 2:
                        assert abs(conc_val - 100) < 0.01
                    else:
                        assert abs(conc_val) < 0.01  # M3 should be 0
        finally:
            output_path.unlink()

    def test_ibot_with_tbnpolys_output(self):
        """Test IBOT output generation."""
        tbn = create_test_tbn()
        polymers = [np.array([1, 0, 0]), np.array([0, 1, 0]), np.array([1, 1, 0])]
        on_target_indices = {0}

        # Create a simple reaction
        reaction_vec = np.array([-1, 1, 0])
        reactions = [Reaction(reaction_vec)]

        ibot = IBOTAlgorithm(tbn, polymers, on_target_indices, reactions)

        # Run the algorithm
        concentration_exponents = ibot.run()

        # Check that only assigned polymers are returned
        # Polymer 0 is on-target (μ=1), polymer 1 gets assigned via reaction
        # Polymer 2 is never involved in any reaction, so it remains unassigned (μ=0)
        assert len(concentration_exponents) == 2
        assert 0 in concentration_exponents  # On-target
        assert 1 in concentration_exponents  # Assigned via reaction
        assert 2 not in concentration_exponents  # Unassigned, removed

        # Test output generation (just check it doesn't crash)
        with tempfile.NamedTemporaryFile(suffix=".tbnpolys", delete=False) as f:
            output_path = Path(f.name)

        try:
            ibot.generate_tbnpolys_output(output_path)
            assert output_path.exists()

            # Read and verify content
            content = output_path.read_text()
            assert "ON-TARGET POLYMERS" in content
            assert "OFF-TARGET POLYMERS" in content
            assert "# μ:" in content
        finally:
            output_path.unlink()


class TestEndToEnd:
    """End-to-end integration tests."""

    def test_and_gate_example(self):
        """Test with the and_gate example files."""
        # Check if example files exist
        tbn_file = Path("extensions/my_inputs/and_gate.tbn")
        on_target_file = Path("extensions/my_inputs/and_gate_on-target.tbnpolys")

        if not tbn_file.exists() or not on_target_file.exists():
            pytest.skip("Example files not found")

        # Parse TBN
        monomers, binding_site_index, concentration_units, _ = TBNParser.parse_file(str(tbn_file))
        assert concentration_units is None  # Should have no concentrations

        tbn = TBN(monomers, binding_site_index, concentration_units)

        # Compute polymer basis
        basis_computer = PolymerBasisComputer(tbn)
        polymers = basis_computer.compute_polymer_basis()
        polymer_vectors = [p.monomer_counts for p in polymers]

        assert len(polymer_vectors) > 0

        # Set up canonical reactions
        reactions_computer = CanonicalReactionsComputer(tbn)
        on_target_indices = reactions_computer.load_on_target_polymers(on_target_file, polymer_vectors)

        assert len(on_target_indices) == 4  # Based on the example file

        reactions_computer.setup_matrices(polymer_vectors, on_target_indices)
        reactions = reactions_computer.compute_irreducible_canonical_reactions()

        assert len(reactions) > 0

        # Check detailed balance
        violating = reactions_computer.check_on_target_detailed_balance(reactions)
        assert violating is None  # Should be balanced

        # Run IBOT
        ibot = IBOTAlgorithm(tbn, polymer_vectors, on_target_indices, reactions)
        concentration_exponents = ibot.run()

        # Check results
        # Only assigned polymers should be in the result
        assert len(concentration_exponents) <= len(polymer_vectors)
        # All on-target polymers should be included and have μ = 1
        for idx in on_target_indices:
            assert idx in concentration_exponents
            assert concentration_exponents[idx] == 1.0  # On-target should have μ = 1
        # Some off-target polymers should have been assigned
        assert len(concentration_exponents) > len(on_target_indices)
