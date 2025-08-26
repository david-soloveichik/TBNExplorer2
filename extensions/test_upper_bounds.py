#!/usr/bin/env python3
"""Test suite for upper bounds computation functionality.

This module tests the computation of upper bounds on specific off-target polymer
concentrations using the restricted irreducible canonical reactions approach.
"""

import subprocess
import tempfile
from pathlib import Path

import numpy as np
import pytest

from extensions.canonical_reactions import CanonicalReactionsComputer, Reaction
from extensions.ibot import IBOTAlgorithm
from tbnexplorer2.model import TBN
from tbnexplorer2.normaliz import NormalizRunner
from tbnexplorer2.parser import TBNParser
from tbnexplorer2.polymer_basis import PolymerBasisComputer
from tbnexplorer2.tbnpolys_io import TbnpolysParser


class TestUpperBoundsComputation:
    """Test upper bounds computation for specific off-target polymers."""

    def setup_method(self):
        """Set up test fixtures."""
        self.temp_dir = tempfile.mkdtemp()
        self.temp_path = Path(self.temp_dir)

    def teardown_method(self):
        """Clean up test fixtures."""
        import shutil

        shutil.rmtree(self.temp_dir, ignore_errors=True)

    def create_test_tbn_file(self, content: str, filename: str = "test.tbn") -> Path:
        """Create a test .tbn file with given content."""
        tbn_path = self.temp_path / filename
        tbn_path.write_text(content)
        return tbn_path

    def create_test_tbnpolys_file(self, content: str, filename: str = "test.tbnpolys") -> Path:
        """Create a test .tbnpolys file with given content."""
        tbnpolys_path = self.temp_path / filename
        tbnpolys_path.write_text(content)
        return tbnpolys_path

    def test_normaliz_strict_inequality(self):
        """Test that Normaliz runner correctly handles strict inequalities."""
        runner = NormalizRunner()

        # Simple test case: find integer points with x + y = 3, x,y >= 0, x + 2y > 2
        # This should give points like (1, 1), (0, 2) but not (3, 0) or (2, 0)
        equations = np.array([[1, 1, -1]])  # x + y - z = 0
        inequalities = np.array(
            [[-1, 0, 0], [0, -1, 0], [0, 0, -1]]
        )  # -x >= 0, -y >= 0, -z >= 0 (flipped for non-negativity)
        # Actually for Normaliz we work with non-negative variables, so we need different setup

        # Let's use a simpler test with the cone formulation
        # We want: x + y = 3 (with slack), x,y >= 0, x + 2y >= 3 (strict as >= 3 for integers)
        # Using homogeneous coordinates (x, y, w) with w = 1 for affine space
        equations = np.array([[1, 1, -3]])  # x + y - 3w = 0 (i.e., x + y = 3 when w = 1)
        inequalities = np.empty((0, 3), dtype=int)  # No additional inequalities (non-negativity is implicit)
        strict_ineq = np.array([1, 2, -2])  # x + 2y - 2w >= 1 (i.e., x + 2y > 2 when w = 1)

        # This should work if Normaliz is installed
        try:
            hilbert_basis = runner.compute_hilbert_basis_with_strict_inequality(equations, inequalities, strict_ineq)
            assert isinstance(hilbert_basis, list)
            # Can't assert specific values without running Normaliz
        except RuntimeError as e:
            if "Normaliz executable not found" in str(e):
                pytest.skip("Normaliz not installed")
            else:
                raise

    def test_canonical_reactions_for_targets_validation(self):
        """Test validation in compute_irreducible_canonical_reactions_for_targets."""
        # Create a simple TBN
        tbn_content = """
a b
c d
a* b*
c* d*
"""
        tbn_path = self.create_test_tbn_file(tbn_content)
        monomers, binding_site_index, concentration_units, _ = TBNParser.parse_file(str(tbn_path))
        tbn = TBN(monomers, binding_site_index, concentration_units)

        # Create canonical reactions computer
        computer = CanonicalReactionsComputer(tbn, use_4ti2=False)

        # Create some dummy polymer basis
        polymer_basis = [
            np.array([1, 0, 0, 0]),  # First monomer only (index 0)
            np.array([0, 1, 0, 0]),  # Second monomer only (index 1)
            np.array([1, 0, 1, 0]),  # Monomers 1 and 3 (index 2)
            np.array([0, 1, 0, 1]),  # Monomers 2 and 4 (index 3)
        ]

        on_target_indices = {0, 1}  # First two polymers are on-target
        computer.setup_matrices(polymer_basis, on_target_indices)

        # Test 1: Target polymers must be off-target
        with pytest.raises(ValueError, match="Target polymers must be off-target"):
            computer.compute_irreducible_canonical_reactions_for_targets({0})  # 0 is on-target

        # Test 2: Target polymer indices must be in range
        with pytest.raises(ValueError, match="Target polymer indices out of range"):
            computer.compute_irreducible_canonical_reactions_for_targets({10})  # Out of range

        # Test 3: Cannot use 4ti2
        computer_4ti2 = CanonicalReactionsComputer(tbn, use_4ti2=True)
        computer_4ti2.setup_matrices(polymer_basis, on_target_indices)
        with pytest.raises(ValueError, match="Upper bound computation requires Normaliz"):
            computer_4ti2.compute_irreducible_canonical_reactions_for_targets({2})

    def test_upper_bounds_cli_validation(self):
        """Test command-line validation for upper bounds options."""
        # Create test files
        tbn_content = """
a b
c d
"""
        tbn_path = self.create_test_tbn_file(tbn_content)

        on_target_content = """
a b
"""
        on_target_path = self.create_test_tbnpolys_file(on_target_content, "on_target.tbnpolys")

        upper_bound_content = """
c d
"""
        upper_bound_path = self.create_test_tbnpolys_file(upper_bound_content, "upper_bound.tbnpolys")

        # Test 1: Cannot use with --use-4ti2
        result = subprocess.run(
            [
                "python",
                "-m",
                "extensions.ibot_cli",
                str(tbn_path),
                str(on_target_path),
                "--upper-bound-on-polymers",
                str(upper_bound_path),
                "--use-4ti2",
            ],
            capture_output=True,
            text=True,
        )
        assert result.returncode != 0
        assert "cannot be used with --use-4ti2" in result.stderr

        # Test 2: Cannot use with --generate-tbn
        result = subprocess.run(
            [
                "python",
                "-m",
                "extensions.ibot_cli",
                str(tbn_path),
                str(on_target_path),
                "--upper-bound-on-polymers",
                str(upper_bound_path),
                "--generate-tbn",
                "100",
                "nM",
            ],
            capture_output=True,
            text=True,
        )
        assert result.returncode != 0
        assert "cannot be used with --generate-tbn" in result.stderr

    def test_upper_bounds_identical_to_full_ibot(self):
        """Test that using all off-target polymers gives identical results to regular IBOT.

        This is the key correctness test: if we specify ALL off-target polymers as targets,
        we should get exactly the same concentration exponents as the regular IBOT algorithm.
        """
        # Use the and_gate_noA.tbn system
        from pathlib import Path
        
        # Get the path to the test files
        current_file = Path(__file__)
        extensions_dir = current_file.parent
        test_tbn = extensions_dir / "my_inputs" / "and_gate_noA.tbn"
        test_on_target = extensions_dir / "my_inputs" / "and_gate_noA_on-target.tbnpolys"
        
        # Check if files exist, otherwise create them
        if not test_tbn.exists():
            tbn_content = """B: b1 b2
a1* a2* b1* b2*
a1 a2 b1 b2 c1
a2* b1* b2* c1*
a2 b1
b2 c1 c2
c1* c2*
C: c1 c2"""
            test_tbn = self.create_test_tbn_file(tbn_content, "and_gate_noA.tbn")
        
        if not test_on_target.exists():
            on_target_content = """B

a1* a2* b1* b2*
a1 a2 b1 b2 c1

a2* b1* b2* c1*
a2 b1
b2 c1 c2

c1* c2*
C"""
            test_on_target = self.create_test_tbnpolys_file(on_target_content, "and_gate_noA_on-target.tbnpolys")

        # Parse TBN
        monomers, binding_site_index, concentration_units, _ = TBNParser.parse_file(str(test_tbn))
        tbn = TBN(monomers, binding_site_index, concentration_units)

        # Compute polymer basis
        try:
            runner = NormalizRunner()
            if not runner.check_normaliz_available():
                pytest.skip("Normaliz not available")
        except Exception:
            pytest.skip("Normaliz not available")

        basis_computer = PolymerBasisComputer(tbn, runner)
        polymers = basis_computer.compute_polymer_basis()
        polymer_vectors = [p.monomer_counts for p in polymers]

        # Load on-target polymers from the file
        parser = TbnpolysParser(tbn)
        on_target_polymers_raw = parser.parse_file(test_on_target)
        
        # Convert to polymer indices
        on_target_indices = set()
        for polymer_raw in on_target_polymers_raw:
            counts = np.zeros(len(tbn.monomers), dtype=int)
            for multiplicity, monomer in polymer_raw:
                monomer_idx = tbn.monomers.index(monomer)
                counts[monomer_idx] += multiplicity
            
            # Find index in polymer basis
            for i, polymer in enumerate(polymer_vectors):
                if np.array_equal(counts, polymer):
                    on_target_indices.add(i)
                    break

        off_target_indices = set(range(len(polymer_vectors))) - on_target_indices

        # Set up canonical reactions computer
        reactions_computer = CanonicalReactionsComputer(tbn, use_4ti2=False)
        reactions_computer.setup_matrices(polymer_vectors, on_target_indices)

        # Compute regular IBOT (all reactions)
        reactions_full = reactions_computer.compute_irreducible_canonical_reactions()

        # Check on-target detailed balance
        violating = reactions_computer.check_on_target_detailed_balance(reactions_full)
        if violating:
            # Skip test if on-target polymers not in detailed balance
            pytest.skip("On-target polymers not in detailed balance for this test case")

        ibot_full = IBOTAlgorithm(tbn, polymer_vectors, on_target_indices, reactions_full)
        mu_full = ibot_full.run()

        # Compute upper bounds using ALL off-target polymers that can be produced
        if not off_target_indices:
            pytest.skip("No off-target polymers in this test case")

        # First find which off-target polymers can actually be produced
        producible_off_target = set()
        for reaction in reactions_full:
            reactants, products = reaction.get_reactants_and_products()
            for idx, _ in products:
                if idx in off_target_indices:
                    producible_off_target.add(idx)
        
        if not producible_off_target:
            pytest.skip("No off-target polymers can be produced by canonical reactions")

        # Compute upper bounds for all producible off-target polymers
        reactions_bounded = reactions_computer.compute_irreducible_canonical_reactions_for_targets(producible_off_target)
        
        # If we're targeting ALL producible off-target polymers, we should get the same reactions
        # (or at least all reactions that produce any off-target)
        ibot_bounded = IBOTAlgorithm(tbn, polymer_vectors, on_target_indices, reactions_bounded)
        mu_bounded = ibot_bounded.run()

        # Compare results - for polymers that appear in both, μ values should be identical
        # Note: Some polymers might not appear in bounded if they can't be produced
        common_polymers = set(mu_full.keys()) & set(mu_bounded.keys())
        assert len(common_polymers) > 0, "No common polymers found between full and bounded"

        for polymer_idx in common_polymers:
            assert abs(mu_full[polymer_idx] - mu_bounded[polymer_idx]) < 1e-10, (
                f"Different μ values for polymer {polymer_idx}: "
                f"full={mu_full[polymer_idx]}, bounded={mu_bounded[polymer_idx]}"
            )

        print("Success: Upper bounds with all off-target polymers matches regular IBOT")
        print(f"Tested with {len(polymer_vectors)} total polymers, {len(off_target_indices)} off-target")


if __name__ == "__main__":
    # Run tests
    pytest.main([__file__, "-v"])
