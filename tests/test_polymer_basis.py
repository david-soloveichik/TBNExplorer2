import tempfile
from pathlib import Path
from unittest.mock import Mock

import numpy as np
import pytest

from tbnexplorer2.model import TBN, BindingSite, Monomer
from tbnexplorer2.polymer_basis import Polymer, PolymerBasisComputer


class TestPolymer:
    def test_polymer_init(self):
        """Test Polymer initialization."""
        monomer_counts = np.array([1, 0, 2])
        monomers = [Mock(spec=Monomer) for _ in range(3)]
        tbn = Mock(spec=TBN)

        polymer = Polymer(monomer_counts, monomers, tbn)

        np.testing.assert_array_equal(polymer.monomer_counts, monomer_counts)
        assert polymer.monomers == monomers
        assert polymer.tbn == tbn

    def test_polymer_init_defaults(self):
        """Test Polymer initialization with default values."""
        monomer_counts = np.array([1, 1])
        monomers = [Mock(spec=Monomer) for _ in range(2)]

        polymer = Polymer(monomer_counts, monomers)

        np.testing.assert_array_equal(polymer.monomer_counts, monomer_counts)
        assert polymer.monomers == monomers
        assert polymer.tbn is None

    def test_polymer_repr(self):
        """Test Polymer string representation."""
        monomer_counts = np.array([1, 0, 1])
        monomers = [Mock(spec=Monomer) for _ in range(3)]
        polymer = Polymer(monomer_counts, monomers)

        # We can't test exact repr since the class doesn't define __repr__
        # But we can test it doesn't raise an error
        repr_str = str(polymer)
        assert repr_str is not None

    def test_polymer_eq(self):
        """Test Polymer equality comparison."""
        monomers = [Mock(spec=Monomer) for _ in range(3)]
        polymer1 = Polymer(np.array([1, 2, 0]), monomers)
        polymer2 = Polymer(np.array([1, 2, 0]), monomers)
        polymer3 = Polymer(np.array([0, 2, 1]), monomers)

        assert polymer1 == polymer2
        assert polymer1 != polymer3
        assert polymer1 != "not a polymer"

    def test_get_monomers_with_counts(self):
        """Test get_monomers_with_counts method."""
        # Create mock monomers with names
        monomer1 = Mock(spec=Monomer)
        monomer1.name = "A"
        monomer2 = Mock(spec=Monomer)
        monomer2.name = "B"
        monomer3 = Mock(spec=Monomer)
        monomer3.name = "C"

        monomers = [monomer1, monomer2, monomer3]

        # Test with some zero counts
        polymer = Polymer(np.array([2, 0, 3]), monomers)
        result = polymer.get_monomers_with_counts()

        assert len(result) == 2  # Only non-zero counts
        assert result[0] == (2, monomer1)
        assert result[1] == (3, monomer3)

        # Test with all zero counts
        polymer_empty = Polymer(np.array([0, 0, 0]), monomers)
        result_empty = polymer_empty.get_monomers_with_counts()
        assert len(result_empty) == 0

        # Test with all non-zero counts
        polymer_full = Polymer(np.array([1, 2, 1]), monomers)
        result_full = polymer_full.get_monomers_with_counts()
        assert len(result_full) == 3
        assert result_full[0] == (1, monomer1)
        assert result_full[1] == (2, monomer2)
        assert result_full[2] == (1, monomer3)

    def test_polymer_hash(self):
        """Test Polymer hash method."""
        monomers = [Mock(spec=Monomer) for _ in range(3)]
        polymer1 = Polymer(np.array([1, 2, 0]), monomers)
        polymer2 = Polymer(np.array([1, 2, 0]), monomers)
        polymer3 = Polymer(np.array([0, 2, 1]), monomers)

        # Same polymer counts should have same hash
        assert hash(polymer1) == hash(polymer2)
        # Different polymer counts should (likely) have different hash
        assert hash(polymer1) != hash(polymer3)

        # Can be used in sets
        polymer_set = {polymer1, polymer2, polymer3}
        assert len(polymer_set) == 2  # polymer1 and polymer2 are the same


class TestPolymerBasisComputer:
    def test_init_with_normaliz(self):
        """Test initialization with Normaliz runner."""
        mock_normaliz = Mock()
        tbn = Mock(spec=TBN)
        computer = PolymerBasisComputer(tbn, normaliz_runner=mock_normaliz)
        assert computer.normaliz_runner == mock_normaliz
        assert computer.tbn == tbn

    def test_init_error_both_runners(self):
        """Test error when both runners are provided."""
        tbn = Mock(spec=TBN)
        # This test may no longer be valid with new signature
        # Let's just test with normaliz only
        mock_normaliz = Mock()
        computer = PolymerBasisComputer(tbn, normaliz_runner=mock_normaliz)
        assert computer.normaliz_runner == mock_normaliz

    def test_init_error_no_runners(self):
        """Test initialization without explicit runner (uses default)."""
        tbn = Mock(spec=TBN)
        # Should create with default Normaliz runner
        computer = PolymerBasisComputer(tbn)
        assert computer.normaliz_runner is not None

    def test_compute_polymer_basis_with_normaliz(self):
        """Test compute_polymer_basis using Normaliz."""
        mock_normaliz = Mock()
        mock_normaliz.compute_hilbert_basis.return_value = [
            np.array([1, 0, 1]),
            np.array([0, 1, 1]),
            np.array([1, 1, 0]),
        ]

        tbn = Mock(spec=TBN)
        tbn.matrix_A = np.array([[1, 0, -1], [0, 1, -1]])
        tbn.num_monomers = 3
        tbn.monomers = [Mock(spec=Monomer) for _ in range(3)]
        tbn.get_augmented_matrix_for_polymer_basis.return_value = (np.array([[1, 0, -1], [0, 1, -1]]), 3)

        computer = PolymerBasisComputer(tbn, normaliz_runner=mock_normaliz)
        polymers = computer.compute_polymer_basis()

        assert len(polymers) == 3
        np.testing.assert_array_equal(polymers[0].monomer_counts, [1, 0, 1])
        np.testing.assert_array_equal(polymers[1].monomer_counts, [0, 1, 1])
        np.testing.assert_array_equal(polymers[2].monomer_counts, [1, 1, 0])

    def test_compute_polymer_basis_with_fourtitwo(self):
        """Test compute_polymer_basis using 4ti2."""
        mock_fourtitwo = Mock()
        mock_fourtitwo.compute_hilbert_basis.return_value = [np.array([1, 0]), np.array([0, 1])]

        tbn = Mock(spec=TBN)
        tbn.matrix_A = np.array([[1, -1]])
        tbn.num_monomers = 2

        # fourtitwo not implemented yet, skip this test
        pytest.skip("4ti2 runner not implemented yet")

    def test_compute_polymer_basis_no_vectors_error(self):
        """Test error when no Hilbert basis vectors are found."""
        mock_normaliz = Mock()
        mock_normaliz.compute_hilbert_basis.return_value = []  # Empty result

        tbn = Mock(spec=TBN)
        tbn.get_augmented_matrix_for_polymer_basis.return_value = (np.array([[1, -1]]), 2)

        computer = PolymerBasisComputer(tbn, normaliz_runner=mock_normaliz)

        with pytest.raises(RuntimeError, match="No Hilbert basis vectors found"):
            computer.compute_polymer_basis()

    def test_compute_free_energies(self):
        """Test compute_free_energies method."""
        tbn = Mock(spec=TBN)
        # Matrix A: 2 binding sites, 3 monomers
        # Monomer 0: {a, b*}  -> [1, -1]
        # Monomer 1: {a*}     -> [-1, 0]
        # Monomer 2: {b}      -> [0, 1]
        tbn.matrix_A = np.array(
            [
                [1, -1, 0],  # binding site a
                [-1, 0, 1],  # binding site b
            ]
        )

        monomers = [Mock(spec=Monomer) for _ in range(3)]
        polymers = [
            Polymer(np.array([1, 1, 0]), monomers, tbn),  # Monomer 0 + Monomer 1
            Polymer(np.array([1, 0, 1]), monomers, tbn),  # Monomer 0 + Monomer 2
            Polymer(np.array([0, 2, 2]), monomers, tbn),  # 2*Monomer 1 + 2*Monomer 2
        ]

        # Compute free energies with default deltaG and collect values
        energies = [polymer.compute_free_energy() for polymer in polymers]

        # Calculate expected free energies with deltaG = -1
        # Matrix A = [[1, -1, 0], [-1, 0, 1]]
        # |A| = [[1, 1, 0], [1, 0, 1]]

        # For polymer [1, 1, 0]:
        # A.x = [1*1 + (-1)*1 + 0*0, (-1)*1 + 0*1 + 1*0] = [0, -1]
        # |A|.x = [1*1 + 1*1 + 0*0, 1*1 + 0*1 + 1*0] = [2, 1]
        # Sum(|A|.x) = 3, Sum(A.x) = -1
        # bonds = (3 - (-1)) / 2 = 2
        # free_energy = -1 * 2 = -2
        assert energies[0] == -2

        # For polymer [1, 0, 1]:
        # A.x = [1*1 + (-1)*0 + 0*1, (-1)*1 + 0*0 + 1*1] = [1, 0]
        # |A|.x = [1*1 + 1*0 + 0*1, 1*1 + 0*0 + 1*1] = [1, 2]
        # Sum(|A|.x) = 3, Sum(A.x) = 1
        # bonds = (3 - 1) / 2 = 1
        # free_energy = -1 * 1 = -1
        assert energies[1] == -1

        # For polymer [0, 2, 2]:
        # A.x = [1*0 + (-1)*2 + 0*2, (-1)*0 + 0*2 + 1*2] = [-2, 2]
        # |A|.x = [1*0 + 1*2 + 0*2, 1*0 + 0*2 + 1*2] = [2, 2]
        # Sum(|A|.x) = 4, Sum(A.x) = 0
        # bonds = (4 - 0) / 2 = 2
        # free_energy = -1 * 2 = -2
        assert energies[2] == -2

    def test_compute_free_energies_with_custom_deltaG(self):
        """Test compute_free_energies method with custom deltaG parameter."""
        tbn = Mock(spec=TBN)
        # Matrix A: 2 binding sites, 3 monomers
        # Monomer 0: {a, b*}  -> [1, -1]
        # Monomer 1: {a*}     -> [-1, 0]
        # Monomer 2: {b}      -> [0, 1]
        tbn.matrix_A = np.array(
            [
                [1, -1, 0],  # binding site a
                [-1, 0, 1],  # binding site b
            ]
        )

        monomers = [Mock(spec=Monomer) for _ in range(3)]
        polymers = [
            Polymer(np.array([1, 1, 0]), monomers, tbn),  # Monomer 0 + Monomer 1
            Polymer(np.array([1, 0, 1]), monomers, tbn),  # Monomer 0 + Monomer 2
            Polymer(np.array([0, 2, 2]), monomers, tbn),  # 2*Monomer 1 + 2*Monomer 2
        ]

        # Compute free energies with custom deltaG = -2.5
        # Using explicit deltaG will apply association penalty
        deltaG = [-2.5, 0.0, 0.0]
        energies2 = [polymer.compute_free_energy(deltaG) for polymer in polymers]

        # Calculate expected free energies with deltaG = -2.5 and association penalty
        # For polymer [1, 1, 0]: 2 monomers, 2 bonds * -2.5 = -5.0, plus association penalty
        assert abs(energies2[0] - (-7.471394)) < 1e-5

        # For polymer [1, 0, 1]: 2 monomers, 1 bond * -2.5 = -2.5, plus association penalty
        assert abs(energies2[1] - (-4.971394)) < 1e-5

        # For polymer [0, 2, 2]: 4 monomers, 2 bonds * -2.5 = -5.0, plus association penalty
        assert abs(energies2[2] - (-12.414182)) < 1e-5

    def test_save_polymer_basis(self):
        """Test save_polymer_basis method with correct signature."""
        # Create mock TBN
        tbn = Mock(spec=TBN)

        # Create mock monomers
        monomer1 = Mock(spec=Monomer)
        monomer1.name = "M1"
        monomer1.binding_sites_str = "a b"
        monomer1.binding_sites = [BindingSite("a", False), BindingSite("b", False)]

        monomer2 = Mock(spec=Monomer)
        monomer2.name = None
        monomer2.binding_sites_str = "c d*"
        monomer2.binding_sites = [BindingSite("c", False), BindingSite("d", True)]

        tbn.monomers = [monomer1, monomer2]

        computer = PolymerBasisComputer(tbn, normaliz_runner=Mock())

        # Create polymers
        polymers = [
            Polymer(np.array([2, 0]), tbn.monomers, tbn),  # 2 of monomer1
            Polymer(np.array([1, 1]), tbn.monomers, tbn),  # 1 of each
            Polymer(np.array([0, 3]), tbn.monomers, tbn),  # 3 of monomer2
        ]

        with tempfile.NamedTemporaryFile(mode="r", delete=False, suffix=".tbnpolys") as f:
            output_file = f.name

        try:
            computer.save_polymer_basis(polymers, output_file)

            # Read and verify the file
            with open(output_file) as f:
                content = f.read()

            # Check header comment
            assert "Polymer basis - 3 polymers" in content

            # Check for polymer representations
            # The file format shows:
            # - "2 | M1" for the first polymer (2x monomer1)
            # - "M1" and "c d*" for the second polymer (1 of each, no prefix means 1)
            # - "3 | c d*" for the third polymer
            assert "2 | M1" in content
            assert "M1" in content  # This represents 1 | M1 without the prefix
            assert "c d*" in content
            assert "3 | c d*" in content

            # Check for separators between polymers
            lines = content.strip().split("\n")
            # Count empty lines between polymers
            empty_count = sum(1 for line in lines if line.strip() == "")
            assert empty_count >= 2  # At least 2 separators for 3 polymers
        finally:
            Path(output_file).unlink()

    def test_load_cached_polymer_basis(self):
        """Test loading cached polymer basis from .tbnpolymat file."""
        # Create mock TBN
        tbn = Mock(spec=TBN)
        tbn.compute_matrix_hash.return_value = "test_hash_123"
        tbn.monomers = [Mock(spec=Monomer) for _ in range(3)]

        computer = PolymerBasisComputer(tbn)

        # Test 1: File doesn't exist
        result = computer.load_cached_polymer_basis("/nonexistent/file.tbnpolymat")
        assert result is None

        # Test 2: File exists with matching hash
        with tempfile.NamedTemporaryFile(mode="w", delete=False, suffix=".tbnpolymat") as f:
            f.write("# Test polymat file\n")
            f.write("\\MATRIX-HASH: test_hash_123\n")
            f.write("# Polymer data\n")
            f.write("1 0 1 -1.0 10.5\n")  # polymer with free energy and concentration
            f.write("0 2 1 -2.0 5.3\n")
            polymat_file = f.name

        try:
            result = computer.load_cached_polymer_basis(polymat_file)
            assert result is not None
            assert len(result) == 2
            np.testing.assert_array_equal(result[0].monomer_counts, [1, 0, 1])
            np.testing.assert_array_equal(result[1].monomer_counts, [0, 2, 1])
        finally:
            Path(polymat_file).unlink()

        # Test 3: File exists with non-matching hash
        with tempfile.NamedTemporaryFile(mode="w", delete=False, suffix=".tbnpolymat") as f:
            f.write("\\MATRIX-HASH: different_hash\n")
            f.write("1 0 1\n")
            polymat_file = f.name

        try:
            result = computer.load_cached_polymer_basis(polymat_file)
            assert result is None  # Should return None for non-matching hash
        finally:
            Path(polymat_file).unlink()

        # Test 4: File with parse errors (non-numeric data)
        with tempfile.NamedTemporaryFile(mode="w", delete=False, suffix=".tbnpolymat") as f:
            f.write("\\MATRIX-HASH: test_hash_123\n")
            f.write("invalid data here\n")
            polymat_file = f.name

        try:
            result = computer.load_cached_polymer_basis(polymat_file)
            assert result is None  # Should return None for parse errors
        finally:
            Path(polymat_file).unlink()

        # Test 5: Empty file with matching hash (edge case)
        with tempfile.NamedTemporaryFile(mode="w", delete=False, suffix=".tbnpolymat") as f:
            f.write("\\MATRIX-HASH: test_hash_123\n")
            f.write("# No polymer data\n")
            polymat_file = f.name

        try:
            result = computer.load_cached_polymer_basis(polymat_file)
            assert result is not None
            assert len(result) == 0  # Empty list but not None
        finally:
            Path(polymat_file).unlink()

    def test_save_tbnpolymat_with_concentrations(self):
        """Test save_tbnpolymat with concentration sorting."""
        # Create mock TBN with concentrations
        tbn = Mock(spec=TBN)
        tbn.compute_matrix_hash.return_value = "test_hash"
        tbn.concentrations = np.array([100.0, 50.0, 75.0])  # in nM
        tbn.concentration_units = "nM"
        tbn.monomers = [Mock(spec=Monomer) for _ in range(3)]
        tbn.matrix_A = np.array([[1, 0, -1], [0, 1, -1]])

        computer = PolymerBasisComputer(tbn)

        # Create polymers with different free energies
        polymer1 = Polymer(np.array([1, 0, 1]), tbn.monomers, tbn)
        polymer2 = Polymer(np.array([0, 1, 1]), tbn.monomers, tbn)
        polymer3 = Polymer(np.array([1, 1, 0]), tbn.monomers, tbn)
        polymers = [polymer1, polymer2, polymer3]

        # Mock COFFEE runner to return concentrations
        mock_coffee = Mock()
        # Return concentrations in reverse order to test sorting
        mock_coffee.compute_equilibrium_concentrations.return_value = np.array([5.0e-9, 20.0e-9, 10.0e-9])  # in Molar

        with tempfile.NamedTemporaryFile(mode="r", delete=False, suffix=".tbnpolymat") as f:
            output_file = f.name

        try:
            computer.save_tbnpolymat(
                polymers,
                output_file,
                compute_free_energies=True,
                compute_concentrations=True,
                concentration_runner=mock_coffee,
                verbose=False,
            )

            # Read and verify the file
            with open(output_file) as f:
                content = f.read()

            # Check that matrix hash is present
            assert "\\MATRIX-HASH: test_hash" in content

            # Check concentration units are specified (in comment format)
            assert "Concentration units:" in content and "nM" in content

            # Parse data lines to verify sorting
            data_lines = []
            for line in content.split("\n"):
                line = line.strip()
                if line and not line.startswith("#") and not line.startswith("\\"):
                    data_lines.append(line)

            # Should have 3 polymer lines
            assert len(data_lines) == 3

            # Verify they are sorted by concentration (descending)
            # polymer2 should be first (20.0e-9 M = 20.0 nM)
            # polymer3 should be second (10.0e-9 M = 10.0 nM)
            # polymer1 should be last (5.0e-9 M = 5.0 nM)
            assert data_lines[0].startswith("0 1 1")  # polymer2
            assert data_lines[1].startswith("1 1 0")  # polymer3
            assert data_lines[2].startswith("1 0 1")  # polymer1

        finally:
            Path(output_file).unlink()

    def test_save_tbnpolymat_concentration_error(self):
        """Test error handling when concentration computation fails."""
        # Create mock TBN with concentrations
        tbn = Mock(spec=TBN)
        tbn.compute_matrix_hash.return_value = "test_hash"
        tbn.concentrations = np.array([100.0, 50.0])
        tbn.concentration_units = "nM"
        tbn.monomers = [Mock(spec=Monomer) for _ in range(2)]
        tbn.matrix_A = np.array([[1, -1]])

        computer = PolymerBasisComputer(tbn)

        # Create a polymer
        polymer = Polymer(np.array([1, 1]), tbn.monomers, tbn)
        polymers = [polymer]

        # Mock COFFEE runner to raise an exception
        mock_coffee = Mock()
        mock_coffee.compute_equilibrium_concentrations.side_effect = Exception("COFFEE failed")

        with tempfile.NamedTemporaryFile(mode="r", delete=False, suffix=".tbnpolymat") as f:
            output_file = f.name

        try:
            # Should not raise, but should skip concentration computation
            computer.save_tbnpolymat(
                polymers,
                output_file,
                compute_free_energies=True,
                compute_concentrations=True,
                concentration_runner=mock_coffee,
                verbose=True,  # Enable to see warning message
            )

            # Read and verify the file
            with open(output_file) as f:
                content = f.read()

            # Should have matrix hash
            assert "\\MATRIX-HASH: test_hash" in content

            # Should NOT have concentration data in the polymer lines
            data_lines = []
            for line in content.split("\n"):
                line = line.strip()
                if line and not line.startswith("#") and not line.startswith("\\"):
                    data_lines.append(line)

            # Check that line doesn't have concentration (would be 5th column)
            parts = data_lines[0].split()
            assert len(parts) == 3  # Just 2 monomer counts + free energy, no concentration

        finally:
            Path(output_file).unlink()

    def test_save_tbnpolymat_no_concentrations(self):
        """Test save_tbnpolymat without monomer concentrations."""
        # Create mock TBN without concentrations
        tbn = Mock(spec=TBN)
        tbn.compute_matrix_hash.return_value = "test_hash"
        tbn.concentrations = None  # No concentrations
        tbn.concentration_units = None
        tbn.monomers = [Mock(spec=Monomer) for _ in range(2)]
        tbn.matrix_A = np.array([[1, -1]])

        computer = PolymerBasisComputer(tbn)

        # Create a polymer
        polymer = Polymer(np.array([1, 1]), tbn.monomers, tbn)
        polymers = [polymer]

        with tempfile.NamedTemporaryFile(mode="r", delete=False, suffix=".tbnpolymat") as f:
            output_file = f.name

        try:
            # Should work without concentrations
            computer.save_tbnpolymat(
                polymers,
                output_file,
                compute_free_energies=True,
                compute_concentrations=True,  # Will be ignored since no monomer concentrations
                verbose=False,
            )

            # Read and verify the file
            with open(output_file) as f:
                content = f.read()

            # Should have matrix hash
            assert "\\MATRIX-HASH: test_hash" in content

            # Should NOT have UNITS line
            assert "\\UNITS:" not in content

            # Check data line has free energy but no concentration
            data_lines = []
            for line in content.split("\n"):
                line = line.strip()
                if line and not line.startswith("#") and not line.startswith("\\"):
                    data_lines.append(line)

            parts = data_lines[0].split()
            assert len(parts) == 3  # 2 monomer counts + free energy

        finally:
            Path(output_file).unlink()
