"""
Extensive tests for matrix hashing functionality to ensure reliability of caching.
"""

import tempfile
from pathlib import Path

import numpy as np

from tbnexplorer2.model import TBN, BindingSite, Monomer
from tbnexplorer2.polymer_basis import PolymerBasisComputer


class TestMatrixHashingExtensive:
    """Extensive tests for matrix hash computation edge cases."""

    def test_hash_stability_across_runs(self):
        """Test that hash is stable across multiple computations."""
        binding_sites = [
            BindingSite("a", False),
            BindingSite("a", True),
            BindingSite("b", False),
            BindingSite("b", True),
            BindingSite("c", False),
            BindingSite("c", True),
        ]

        monomers = [
            Monomer("M1", [binding_sites[0], binding_sites[1], binding_sites[2]], None, "a a* b"),
            Monomer("M2", [binding_sites[3], binding_sites[4], binding_sites[5]], None, "b* c c*"),
            Monomer("M3", [binding_sites[0], binding_sites[3], binding_sites[4]], None, "a b* c"),
        ]

        binding_site_index = {"a": 0, "a*": 1, "b": 2, "b*": 3, "c": 4, "c*": 5}

        tbn = TBN(monomers, binding_site_index)

        # Compute hash multiple times
        hashes = [tbn.compute_matrix_hash() for _ in range(10)]

        # All hashes should be identical
        assert len(set(hashes)) == 1
        assert all(h == hashes[0] for h in hashes)

    def test_hash_different_for_transposed_matrix(self):
        """Test that transposed matrix produces different hash."""
        # Create a non-square matrix
        binding_site_index = {"a": 0, "b": 1, "c": 2}

        monomers1 = [
            Monomer("M1", [BindingSite("a", False)], None, "a"),
            Monomer("M2", [BindingSite("b", False)], None, "b"),
        ]

        tbn1 = TBN(monomers1, binding_site_index)
        hash1 = tbn1.compute_matrix_hash()

        # Create TBN with different shape (more monomers)
        monomers2 = [
            Monomer("M1", [BindingSite("a", False)], None, "a"),
            Monomer("M2", [BindingSite("b", False)], None, "b"),
            Monomer("M3", [BindingSite("c", False)], None, "c"),
        ]

        tbn2 = TBN(monomers2, binding_site_index)
        hash2 = tbn2.compute_matrix_hash()

        assert hash1 != hash2

    def test_hash_sensitive_to_sign_changes(self):
        """Test that hash changes when signs in matrix change."""
        binding_site_index = {"a": 0, "a*": 1}

        # TBN with unstar binding site
        monomers1 = [Monomer("M1", [BindingSite("a", False)], None, "a")]
        tbn1 = TBN(monomers1, binding_site_index)
        hash1 = tbn1.compute_matrix_hash()

        # TBN with star binding site (negative in matrix)
        monomers2 = [Monomer("M1", [BindingSite("a", True)], None, "a*")]
        tbn2 = TBN(monomers2, binding_site_index)
        hash2 = tbn2.compute_matrix_hash()

        assert hash1 != hash2

    def test_hash_sensitive_to_zero_columns(self):
        """Test that hash is different when zero columns are added."""
        binding_site_index = {"a": 0, "b": 1}

        # TBN with one monomer
        monomers1 = [Monomer("M1", [BindingSite("a", False)], None, "a")]
        tbn1 = TBN(monomers1, binding_site_index)
        hash1 = tbn1.compute_matrix_hash()

        # TBN with two monomers (second has all zeros)
        monomers2 = [
            Monomer("M1", [BindingSite("a", False)], None, "a"),
            Monomer("M2", [], None, ""),  # Empty monomer creates zero column
        ]
        tbn2 = TBN(monomers2, binding_site_index)
        hash2 = tbn2.compute_matrix_hash()

        assert hash1 != hash2

    def test_hash_order_sensitivity(self):
        """Test that monomer order affects hash (columns are ordered)."""
        binding_site_index = {"a": 0, "b": 1}

        # First order
        monomers1 = [
            Monomer("M1", [BindingSite("a", False)], None, "a"),
            Monomer("M2", [BindingSite("b", False)], None, "b"),
        ]
        tbn1 = TBN(monomers1, binding_site_index)
        hash1 = tbn1.compute_matrix_hash()

        # Reversed order
        monomers2 = [
            Monomer("M2", [BindingSite("b", False)], None, "b"),
            Monomer("M1", [BindingSite("a", False)], None, "a"),
        ]
        tbn2 = TBN(monomers2, binding_site_index)
        hash2 = tbn2.compute_matrix_hash()

        assert hash1 != hash2

    def test_hash_with_large_matrix(self):
        """Test hash computation with a large matrix."""
        # Create a large TBN with many binding sites and monomers
        n_sites = 50
        n_monomers = 30

        binding_site_index = {}
        for i in range(n_sites):
            binding_site_index[f"s{i}"] = i
            binding_site_index[f"s{i}*"] = i + n_sites

        monomers = []
        for i in range(n_monomers):
            # Create random binding sites for each monomer
            sites = []
            if i % 3 == 0:
                sites.append(BindingSite(f"s{i % n_sites}", False))
            if i % 3 == 1:
                sites.append(BindingSite(f"s{i % n_sites}", True))
            if i % 3 == 2:
                sites.extend([BindingSite(f"s{i % n_sites}", False), BindingSite(f"s{(i + 1) % n_sites}", True)])

            monomers.append(Monomer(f"M{i}", sites, None, f"monomer_{i}"))

        tbn = TBN(monomers, binding_site_index)

        # Should compute hash without errors
        hash1 = tbn.compute_matrix_hash()
        hash2 = tbn.compute_matrix_hash()

        assert hash1 == hash2
        assert len(hash1) == 64  # SHA256 always produces 64 hex chars

    def test_hash_with_empty_matrix(self):
        """Test hash computation with empty matrix (no monomers)."""
        binding_site_index = {"a": 0, "b": 1}

        # TBN with no monomers
        tbn = TBN([], binding_site_index)

        # Should still compute a valid hash
        hash_val = tbn.compute_matrix_hash()
        assert len(hash_val) == 64

        # Hash should be consistent
        assert hash_val == tbn.compute_matrix_hash()

    def test_hash_numeric_precision(self):
        """Test that hash uses integer matrix values (no floating point issues)."""
        binding_site_index = {"a": 0, "b": 1}

        # Create two identical TBNs
        monomers1 = [
            Monomer("M1", [BindingSite("a", False)], None, "a"),
            Monomer("M2", [BindingSite("b", False)], None, "b"),
        ]

        monomers2 = [
            Monomer("M1", [BindingSite("a", False)], None, "a"),
            Monomer("M2", [BindingSite("b", False)], None, "b"),
        ]

        tbn1 = TBN(monomers1, binding_site_index)
        tbn2 = TBN(monomers2, binding_site_index)

        # Verify matrices are integer type
        assert tbn1.matrix_A.dtype == np.dtype("int64") or tbn1.matrix_A.dtype == np.dtype("int32")
        assert tbn2.matrix_A.dtype == np.dtype("int64") or tbn2.matrix_A.dtype == np.dtype("int32")

        # Hashes should be identical
        assert tbn1.compute_matrix_hash() == tbn2.compute_matrix_hash()

    def test_hash_with_complex_binding_patterns(self):
        """Test hash with complex binding patterns including self-binding."""
        binding_site_index = {"a": 0, "b": 1, "c": 2}

        # Complex monomers with various binding patterns
        monomers = [
            # Self-binding monomer
            Monomer("Self", [BindingSite("a", False), BindingSite("a", True)], None, "a a*"),
            # Multiple same sites
            Monomer(
                "Multi", [BindingSite("b", False), BindingSite("b", False), BindingSite("b", True)], None, "b b b*"
            ),
            # Mixed sites
            Monomer(
                "Mixed",
                [BindingSite("a", False), BindingSite("b", True), BindingSite("c", False), BindingSite("c", True)],
                None,
                "a b* c c*",
            ),
        ]

        tbn = TBN(monomers, binding_site_index)

        # Check matrix values are as expected
        matrix = tbn.matrix_A

        # Self-binding monomer should cancel (net 0 for site 'a')
        assert matrix[0, 0] == 0  # a and a* cancel out (1 + (-1) = 0)

        # Multi monomer should have net +1 for b
        assert matrix[1, 1] == 1  # Two b sites minus one b* site (2 + (-1) = 1)

        # Hash should be consistent
        hash1 = tbn.compute_matrix_hash()
        hash2 = tbn.compute_matrix_hash()
        assert hash1 == hash2


class TestCachingWorkflowExtensive:
    """Test the complete caching workflow with various scenarios."""

    def create_test_tbn(self, seed=42):
        """Create a test TBN with deterministic randomness."""
        np.random.seed(seed)

        n_sites = 6
        binding_site_index = {}
        binding_sites = []

        for i in range(n_sites):
            name = f"s{i}"
            binding_site_index[name] = i * 2
            binding_site_index[f"{name}*"] = i * 2 + 1
            binding_sites.append(BindingSite(name, False))
            binding_sites.append(BindingSite(name, True))

        # Create monomers with random binding sites
        monomers = []
        for i in range(4):
            n_bonds = np.random.randint(1, 5)
            sites = []
            for _ in range(n_bonds):
                sites.append(binding_sites[np.random.randint(0, len(binding_sites))])

            sites_str = " ".join(str(s) for s in sites)
            monomers.append(Monomer(f"M{i}", sites, 10.0 * (i + 1), f"{sites_str}, {10.0 * (i + 1)}"))

        return TBN(monomers, binding_site_index, "nM")

    def test_caching_workflow_with_matching_hash(self):
        """Test complete workflow when cached hash matches."""
        tbn = self.create_test_tbn()
        computer = PolymerBasisComputer(tbn, None)

        with tempfile.TemporaryDirectory() as temp_dir:
            polymat_file = Path(temp_dir) / "test.tbnpolymat"

            # Create initial polymers and save
            from tbnexplorer2.polymer_basis import Polymer

            polymers = [
                Polymer(np.array([1, 0, 0, 0]), tbn.monomers, tbn),
                Polymer(np.array([0, 1, 0, 0]), tbn.monomers, tbn),
                Polymer(np.array([1, 1, 0, 0]), tbn.monomers, tbn),
            ]

            computer.save_tbnpolymat(
                polymers, str(polymat_file), compute_free_energies=True, compute_concentrations=False
            )

            # Verify file contains hash
            with open(polymat_file) as f:
                content = f.read()

            expected_hash = tbn.compute_matrix_hash()
            assert f"\\MATRIX-HASH: {expected_hash}" in content

            # Load cached polymers
            cached = computer.load_cached_polymer_basis(str(polymat_file))

            assert cached is not None
            assert len(cached) == 3
            assert np.array_equal(cached[0].monomer_counts, np.array([1, 0, 0, 0]))
            assert np.array_equal(cached[1].monomer_counts, np.array([0, 1, 0, 0]))
            assert np.array_equal(cached[2].monomer_counts, np.array([1, 1, 0, 0]))

    def test_caching_workflow_with_modified_tbn(self):
        """Test that cache is invalidated when TBN changes."""
        tbn1 = self.create_test_tbn(seed=42)
        computer1 = PolymerBasisComputer(tbn1, None)

        with tempfile.TemporaryDirectory() as temp_dir:
            polymat_file = Path(temp_dir) / "test.tbnpolymat"

            # Save with first TBN
            from tbnexplorer2.polymer_basis import Polymer

            polymers = [Polymer(np.array([1, 0, 0, 0]), tbn1.monomers, tbn1)]

            computer1.save_tbnpolymat(
                polymers, str(polymat_file), compute_free_energies=False, compute_concentrations=False
            )

            # Create modified TBN (different seed = different monomers)
            tbn2 = self.create_test_tbn(seed=123)
            computer2 = PolymerBasisComputer(tbn2, None)

            # Hash should be different
            assert tbn1.compute_matrix_hash() != tbn2.compute_matrix_hash()

            # Loading should return None (hash mismatch)
            cached = computer2.load_cached_polymer_basis(str(polymat_file))
            assert cached is None

    def test_caching_with_concentration_changes(self):
        """Test that cache works when only concentrations change (not matrix)."""
        # Create TBN with concentrations
        binding_site_index = {"a": 0, "b": 1}

        monomers1 = [
            Monomer("M1", [BindingSite("a", False)], 100.0, "a, 100"),
            Monomer("M2", [BindingSite("b", False)], 50.0, "b, 50"),
        ]
        tbn1 = TBN(monomers1, binding_site_index, "nM")

        # Create TBN with same structure but different concentrations
        monomers2 = [
            Monomer("M1", [BindingSite("a", False)], 200.0, "a, 200"),
            Monomer("M2", [BindingSite("b", False)], 75.0, "b, 75"),
        ]
        tbn2 = TBN(monomers2, binding_site_index, "nM")

        # Matrix hash should be the same (structure unchanged)
        assert tbn1.compute_matrix_hash() == tbn2.compute_matrix_hash()

        # This allows reusing cached polymer basis even with different concentrations
        computer1 = PolymerBasisComputer(tbn1, None)
        computer2 = PolymerBasisComputer(tbn2, None)

        with tempfile.TemporaryDirectory() as temp_dir:
            polymat_file = Path(temp_dir) / "test.tbnpolymat"

            # Save with first TBN
            from tbnexplorer2.polymer_basis import Polymer

            polymers = [Polymer(np.array([1, 0]), tbn1.monomers, tbn1)]

            computer1.save_tbnpolymat(
                polymers, str(polymat_file), compute_free_energies=False, compute_concentrations=False
            )

            # Load with second TBN (different concentrations, same structure)
            cached = computer2.load_cached_polymer_basis(str(polymat_file))

            assert cached is not None
            assert len(cached) == 1

    def test_caching_robustness_to_file_corruption(self):
        """Test that caching gracefully handles corrupted files."""
        tbn = self.create_test_tbn()
        computer = PolymerBasisComputer(tbn, None)

        with tempfile.TemporaryDirectory() as temp_dir:
            polymat_file = Path(temp_dir) / "test.tbnpolymat"

            # Test various corruption scenarios

            # 1. Malformed hash line
            with open(polymat_file, "w") as f:
                f.write("\\MATRIX-HASH: not_a_valid_hash\n")
                f.write("1 0 0 0\n")

            cached = computer.load_cached_polymer_basis(str(polymat_file))
            assert cached is None  # Should return None, not crash

            # 2. Invalid polymer data
            with open(polymat_file, "w") as f:
                f.write(f"\\MATRIX-HASH: {tbn.compute_matrix_hash()}\n")
                f.write("not_a_number 0 0 0\n")

            cached = computer.load_cached_polymer_basis(str(polymat_file))
            assert cached is None  # Should return None, not crash

            # 3. Wrong number of monomer counts
            with open(polymat_file, "w") as f:
                f.write(f"\\MATRIX-HASH: {tbn.compute_matrix_hash()}\n")
                f.write("1 0\n")  # Only 2 values but TBN has 4 monomers

            cached = computer.load_cached_polymer_basis(str(polymat_file))
            # Should still load (ignores invalid lines)
            assert cached is not None
            assert len(cached) == 0  # No valid polymers loaded

    def test_hash_format_validation(self):
        """Test that hash format is properly validated."""
        tbn = self.create_test_tbn()
        real_hash = tbn.compute_matrix_hash()

        # Valid SHA256 hash should be 64 hex characters
        assert len(real_hash) == 64
        assert all(c in "0123456789abcdef" for c in real_hash)

        computer = PolymerBasisComputer(tbn, None)

        with tempfile.TemporaryDirectory() as temp_dir:
            polymat_file = Path(temp_dir) / "test.tbnpolymat"

            # Test with truncated hash
            with open(polymat_file, "w") as f:
                f.write(f"\\MATRIX-HASH: {real_hash[:32]}\n")  # Only half the hash
                f.write("1 0 0 0\n")

            cached = computer.load_cached_polymer_basis(str(polymat_file))
            assert cached is None  # Should reject truncated hash
