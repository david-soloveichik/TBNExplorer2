import tempfile
from pathlib import Path
from unittest.mock import Mock

import numpy as np
import pytest

from tbnexplorer2.model import TBN, Monomer
from tbnexplorer2.polymer_basis import Polymer, PolymerBasisComputer


class TestPolymer:
    def test_polymer_init(self):
        """Test Polymer initialization."""
        monomer_counts = np.array([1, 0, 2])
        free_energy = -3.5
        concentration = 1.2e-7

        polymer = Polymer(monomer_counts, free_energy, concentration)

        np.testing.assert_array_equal(polymer.monomer_counts, monomer_counts)
        assert polymer.free_energy == free_energy
        assert polymer.concentration == concentration

    def test_polymer_init_defaults(self):
        """Test Polymer initialization with default values."""
        monomer_counts = np.array([1, 1])

        polymer = Polymer(monomer_counts)

        np.testing.assert_array_equal(polymer.monomer_counts, monomer_counts)
        assert polymer.free_energy is None
        assert polymer.concentration is None

    def test_polymer_repr(self):
        """Test Polymer string representation."""
        polymer = Polymer(np.array([1, 0, 1]), -2.0, 1.5e-8)
        repr_str = repr(polymer)

        assert "Polymer" in repr_str
        assert "[1 0 1]" in repr_str
        assert "-2.0" in repr_str
        assert "1.5e-08" in repr_str

    def test_polymer_eq(self):
        """Test Polymer equality comparison."""
        polymer1 = Polymer(np.array([1, 2, 0]))
        polymer2 = Polymer(np.array([1, 2, 0]))
        polymer3 = Polymer(np.array([0, 2, 1]))

        assert polymer1 == polymer2
        assert polymer1 != polymer3
        assert polymer1 != "not a polymer"


class TestPolymerBasisComputer:
    def test_init_with_normaliz(self):
        """Test initialization with Normaliz runner."""
        mock_normaliz = Mock()
        computer = PolymerBasisComputer(normaliz_runner=mock_normaliz)
        assert computer.normaliz_runner == mock_normaliz
        assert computer.fourtitwo_runner is None

    def test_init_with_fourtitwo(self):
        """Test initialization with 4ti2 runner."""
        mock_fourtitwo = Mock()
        computer = PolymerBasisComputer(fourtitwo_runner=mock_fourtitwo)
        assert computer.fourtitwo_runner == mock_fourtitwo
        assert computer.normaliz_runner is None

    def test_init_error_both_runners(self):
        """Test error when both runners are provided."""
        with pytest.raises(ValueError, match="Cannot specify both"):
            PolymerBasisComputer(normaliz_runner=Mock(), fourtitwo_runner=Mock())

    def test_init_error_no_runners(self):
        """Test error when no runners are provided."""
        with pytest.raises(ValueError, match="Must specify either"):
            PolymerBasisComputer()

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

        computer = PolymerBasisComputer(normaliz_runner=mock_normaliz)
        polymers = computer.compute_polymer_basis(tbn)

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

        computer = PolymerBasisComputer(fourtitwo_runner=mock_fourtitwo)
        polymers = computer.compute_polymer_basis(tbn)

        assert len(polymers) == 2
        np.testing.assert_array_equal(polymers[0].monomer_counts, [1, 0])
        np.testing.assert_array_equal(polymers[1].monomer_counts, [0, 1])

    def test_augment_matrix_for_hilbert_basis(self):
        """Test _augment_matrix_for_hilbert_basis method."""
        computer = PolymerBasisComputer(normaliz_runner=Mock())

        tbn = Mock(spec=TBN)
        # Matrix with 2 binding sites, 3 monomers
        tbn.matrix_A = np.array(
            [
                [1, 0, -1],  # binding site 0
                [0, 1, -1],  # binding site 1
            ]
        )

        augmented = computer._augment_matrix_for_hilbert_basis(tbn)

        # Should add singleton star monomers for each binding site
        assert augmented.shape == (2, 5)  # 3 original + 2 singleton
        np.testing.assert_array_equal(augmented[:, 0:3], tbn.matrix_A)
        np.testing.assert_array_equal(augmented[:, 3], [-1, 0])  # singleton for site 0
        np.testing.assert_array_equal(augmented[:, 4], [0, -1])  # singleton for site 1

    def test_remove_fake_monomers(self):
        """Test _remove_fake_monomers method."""
        computer = PolymerBasisComputer(normaliz_runner=Mock())

        # Vectors with 5 entries (3 real monomers, 2 fake)
        vectors = [
            np.array([1, 0, 1, 0, 0]),  # only real monomers
            np.array([0, 1, 0, 1, 0]),  # has fake monomer
            np.array([1, 0, 1, 0, 0]),  # duplicate of first
            np.array([1, 1, 0, 0, 0]),  # only real monomers
        ]

        result = computer._remove_fake_monomers(vectors, 3)

        # Should keep only unique vectors with first 3 entries
        assert len(result) == 3
        np.testing.assert_array_equal(result[0], [1, 0, 1])
        np.testing.assert_array_equal(result[1], [0, 1, 0])
        np.testing.assert_array_equal(result[2], [1, 1, 0])

    def test_compute_free_energies(self):
        """Test compute_free_energies method."""
        computer = PolymerBasisComputer(normaliz_runner=Mock())

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

        polymers = [
            Polymer(np.array([1, 1, 0])),  # Monomer 0 + Monomer 1
            Polymer(np.array([1, 0, 1])),  # Monomer 0 + Monomer 2
            Polymer(np.array([0, 2, 2])),  # 2*Monomer 1 + 2*Monomer 2
        ]

        computer.compute_free_energies(polymers, tbn)

        # Calculate expected free energies
        # For polymer [1, 1, 0]: A.x = [0, -1], |A|.x = [2, 1]
        # bonds = (3 - 1) / 2 = 1
        assert polymers[0].free_energy == -1

        # For polymer [1, 0, 1]: A.x = [1, 0], |A|.x = [1, 2]
        # bonds = (3 - 1) / 2 = 1
        assert polymers[1].free_energy == -1

        # For polymer [0, 2, 2]: A.x = [-2, 2], |A|.x = [2, 2]
        # bonds = (4 - 0) / 2 = 2
        assert polymers[2].free_energy == -2

    def test_save_polymer_basis(self):
        """Test save_polymer_basis method."""
        computer = PolymerBasisComputer(normaliz_runner=Mock())

        # Create mock monomers
        monomer1 = Mock(spec=Monomer)
        monomer1.name = "M1"
        monomer1.binding_sites_str = "a b"

        monomer2 = Mock(spec=Monomer)
        monomer2.name = None
        monomer2.binding_sites_str = "c d*"

        monomers = [monomer1, monomer2]

        # Create polymers
        polymers = [
            Polymer(np.array([2, 0])),  # 2 of monomer1
            Polymer(np.array([1, 1])),  # 1 of each
            Polymer(np.array([0, 3])),  # 3 of monomer2
        ]

        with tempfile.NamedTemporaryFile(mode="r", delete=False, suffix=".txt") as f:
            output_file = f.name

        try:
            computer.save_polymer_basis(polymers, monomers, output_file)

            # Read and verify the file
            with open(output_file) as f:
                content = f.read()

            # Check for polymer representations
            assert "2 | M1" in content
            assert "1 | M1" in content
            assert "1 | c d*" in content
            assert "3 | c d*" in content

            # Check for separators between polymers
            lines = content.strip().split("\n")
            empty_lines = [i for i, line in enumerate(lines) if line.strip() == ""]
            assert len(empty_lines) >= 2  # At least 2 separators for 3 polymers
        finally:
            Path(output_file).unlink()

    def test_compute_polymer_basis_with_output_dir(self):
        """Test compute_polymer_basis with output directory specified."""
        mock_normaliz = Mock()
        mock_normaliz.compute_hilbert_basis.return_value = [np.array([1, 0])]

        tbn = Mock(spec=TBN)
        tbn.matrix_A = np.array([[1, -1]])
        tbn.num_monomers = 2

        computer = PolymerBasisComputer(normaliz_runner=mock_normaliz)

        with tempfile.TemporaryDirectory() as temp_dir:
            computer.compute_polymer_basis(tbn, output_dir=temp_dir)

            # Verify normaliz was called with output_dir
            mock_normaliz.compute_hilbert_basis.assert_called_once()
            call_args = mock_normaliz.compute_hilbert_basis.call_args
            assert call_args[1].get("output_dir") == temp_dir
