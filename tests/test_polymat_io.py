"""
Tests for the polymat_io module.
"""

import os
import tempfile

import numpy as np
import pytest

from tbnexplorer2.polymat_io import (
    PolymatData,
    PolymatReader,
    PolymatWriter,
    check_matrix_hash,
    load_polymat_file,
    save_polymat_file,
)


class TestPolymatData:
    """Tests for PolymatData dataclass."""

    def test_polymat_data_creation(self):
        """Test creating a PolymatData object."""
        polymers = [np.array([1, 0, 2]), np.array([0, 3, 1]), np.array([2, 1, 0])]

        data = PolymatData(
            polymers=polymers,
            n_monomers=3,
            n_polymers=3,
            matrix_hash="test_hash",
            free_energies=np.array([-2.0, -3.5, -1.5]),
            concentrations=np.array([100.0, 50.0, 25.0]),
            concentration_units="nM",
            has_free_energies=True,
            has_concentrations=True,
        )

        assert data.n_monomers == 3
        assert data.n_polymers == 3
        assert data.matrix_hash == "test_hash"
        assert data.concentration_units == "nM"
        assert data.has_free_energies
        assert data.has_concentrations

    def test_get_polymer_data(self):
        """Test getting data for a specific polymer."""
        polymers = [np.array([1, 0, 2]), np.array([0, 3, 1])]

        data = PolymatData(
            polymers=polymers,
            n_monomers=3,
            n_polymers=2,
            free_energies=np.array([-2.0, -3.5]),
            concentrations=np.array([100.0, 50.0]),
            has_free_energies=True,
            has_concentrations=True,
        )

        # Test first polymer
        counts, fe, conc = data.get_polymer_data(0)
        assert np.array_equal(counts, np.array([1, 0, 2]))
        assert fe == -2.0
        assert conc == 100.0

        # Test second polymer
        counts, fe, conc = data.get_polymer_data(1)
        assert np.array_equal(counts, np.array([0, 3, 1]))
        assert fe == -3.5
        assert conc == 50.0

    def test_get_polymer_data_without_optional_fields(self):
        """Test getting polymer data when optional fields are None."""
        polymers = [np.array([1, 0, 2])]

        data = PolymatData(
            polymers=polymers, n_monomers=3, n_polymers=1, has_free_energies=False, has_concentrations=False
        )

        counts, fe, conc = data.get_polymer_data(0)
        assert np.array_equal(counts, np.array([1, 0, 2]))
        assert fe is None
        assert conc is None

    def test_get_polymer_data_out_of_range(self):
        """Test error handling for out of range polymer index."""
        polymers = [np.array([1, 0, 2])]

        data = PolymatData(polymers=polymers, n_monomers=3, n_polymers=1)

        with pytest.raises(IndexError):
            data.get_polymer_data(1)

        with pytest.raises(IndexError):
            data.get_polymer_data(-2)


class TestPolymatWriter:
    """Tests for PolymatWriter class."""

    def test_write_basic_polymat_file(self):
        """Test writing a basic .tbnpolymat file."""
        with tempfile.TemporaryDirectory() as tmpdir:
            file_path = os.path.join(tmpdir, "test.tbnpolymat")

            polymers = [np.array([1, 0, 2]), np.array([0, 3, 1]), np.array([2, 1, 0])]

            data = PolymatData(
                polymers=polymers,
                n_monomers=3,
                n_polymers=3,
                matrix_hash="abc123",
                has_free_energies=False,
                has_concentrations=False,
            )

            writer = PolymatWriter(file_path)
            writer.write(data)

            # Read and verify the file
            with open(file_path) as f:
                lines = f.readlines()

            # Check header
            assert "# TBN Polymer Matrix" in lines[0]
            assert "# Number of polymers: 3" in lines[1]
            assert "# Number of monomers: 3" in lines[2]
            assert "\\MATRIX-HASH: abc123" in lines[3]

            # Check data lines (exclude comments and keyword lines)
            data_lines = [
                line.strip()
                for line in lines
                if not line.startswith("#") and not line.startswith("\\") and line.strip()
            ]
            assert len(data_lines) == 3
            assert data_lines[0] == "1 0 2"
            assert data_lines[1] == "0 3 1"
            assert data_lines[2] == "2 1 0"

    def test_write_with_free_energies(self):
        """Test writing with free energies."""
        with tempfile.TemporaryDirectory() as tmpdir:
            file_path = os.path.join(tmpdir, "test.tbnpolymat")

            polymers = [np.array([1, 0]), np.array([0, 2])]

            data = PolymatData(
                polymers=polymers,
                n_monomers=2,
                n_polymers=2,
                free_energies=np.array([-2.5, -3.0]),
                has_free_energies=True,
                has_concentrations=False,
            )

            writer = PolymatWriter(file_path)
            writer.write(data)

            # Read and verify
            with open(file_path) as f:
                lines = f.readlines()

            # Check columns header
            columns_line = next(line for line in lines if "# Columns:" in line)
            assert "free_energy" in columns_line
            assert "concentration" not in columns_line

            # Check data
            data_lines = [line.strip() for line in lines if not line.startswith("#") and line.strip()]
            assert data_lines[0] == "1 0 -2.5"
            assert data_lines[1] == "0 2 -3.0"

    def test_write_with_concentrations(self):
        """Test writing with concentrations."""
        with tempfile.TemporaryDirectory() as tmpdir:
            file_path = os.path.join(tmpdir, "test.tbnpolymat")

            polymers = [np.array([1, 0]), np.array([0, 2])]

            data = PolymatData(
                polymers=polymers,
                n_monomers=2,
                n_polymers=2,
                free_energies=np.array([-2.5, -3.0]),
                concentrations=np.array([1.5e-7, 2.3e-8]),
                concentration_units="nM",
                has_free_energies=True,
                has_concentrations=True,
            )

            writer = PolymatWriter(file_path)
            writer.write(data)

            # Read and verify
            with open(file_path) as f:
                lines = f.readlines()

            # Check headers
            assert any("# Concentration units: nM" in line for line in lines)

            columns_line = next(line for line in lines if "# Columns:" in line)
            assert "free_energy" in columns_line
            assert "concentration" in columns_line

            # Check data
            data_lines = [line.strip() for line in lines if not line.startswith("#") and line.strip()]
            assert data_lines[0] == "1 0 -2.5 1.50e-07"
            assert data_lines[1] == "0 2 -3.0 2.30e-08"


class TestPolymatReader:
    """Tests for PolymatReader class."""

    def create_test_file(self, tmpdir, content):
        """Helper to create a test .tbnpolymat file."""
        file_path = os.path.join(tmpdir, "test.tbnpolymat")
        with open(file_path, "w") as f:
            f.write(content)
        return file_path

    def test_read_basic_file(self):
        """Test reading a basic .tbnpolymat file."""
        content = """# TBN Polymer Matrix
# Number of polymers: 3
# Number of monomers: 3
\\MATRIX-HASH: test123
# Columns: monomer_counts[1..3]
#
1 0 2
0 3 1
2 1 0
"""

        with tempfile.TemporaryDirectory() as tmpdir:
            file_path = self.create_test_file(tmpdir, content)

            reader = PolymatReader(file_path)
            data = reader.read()

            assert data.n_monomers == 3
            assert data.n_polymers == 3
            assert data.matrix_hash == "test123"
            assert not data.has_free_energies
            assert not data.has_concentrations

            assert len(data.polymers) == 3
            assert np.array_equal(data.polymers[0], np.array([1, 0, 2]))
            assert np.array_equal(data.polymers[1], np.array([0, 3, 1]))
            assert np.array_equal(data.polymers[2], np.array([2, 1, 0]))

    def test_read_with_free_energies(self):
        """Test reading file with free energies."""
        content = """# TBN Polymer Matrix
# Number of polymers: 2
# Number of monomers: 2
# Columns: monomer_counts[1..2] free_energy
#
1 0 -2.5
0 2 -3.0
"""

        with tempfile.TemporaryDirectory() as tmpdir:
            file_path = self.create_test_file(tmpdir, content)

            reader = PolymatReader(file_path)
            data = reader.read()

            assert data.has_free_energies
            assert not data.has_concentrations
            assert data.free_energies is not None
            assert np.allclose(data.free_energies, np.array([-2.5, -3.0]))

    def test_read_with_concentrations(self):
        """Test reading file with free energies and concentrations."""
        content = """# TBN Polymer Matrix
# Number of polymers: 2
# Number of monomers: 2
# Concentration units: nM
# Columns: monomer_counts[1..2] free_energy concentration
#
1 0 -2.5 1.5e-7
0 2 -3.0 2.3e-8
"""

        with tempfile.TemporaryDirectory() as tmpdir:
            file_path = self.create_test_file(tmpdir, content)

            reader = PolymatReader(file_path)
            data = reader.read()

            assert data.has_free_energies
            assert data.has_concentrations
            assert data.concentration_units == "nM"
            assert np.allclose(data.free_energies, np.array([-2.5, -3.0]))
            assert np.allclose(data.concentrations, np.array([1.5e-7, 2.3e-8]))

    def test_read_header_only(self):
        """Test reading only header information."""
        content = """# TBN Polymer Matrix
# Number of polymers: 100
# Number of monomers: 5
\\MATRIX-HASH: xyz789
# Columns: monomer_counts[1..5] free_energy
#
1 0 0 0 0 -1.0
"""

        with tempfile.TemporaryDirectory() as tmpdir:
            file_path = self.create_test_file(tmpdir, content)

            reader = PolymatReader(file_path)
            header = reader.read_header_only()

            assert header["n_monomers"] == 5
            assert header["n_polymers"] == 100
            assert header["matrix_hash"] == "xyz789"
            assert header["has_free_energies"]
            assert not header["has_concentrations"]

    def test_iterate_polymers(self):
        """Test iterating through polymers without loading all into memory."""
        content = """# TBN Polymer Matrix
# Number of polymers: 3
# Number of monomers: 2
# Columns: monomer_counts[1..2] free_energy
#
1 0 -1.0
0 1 -2.0
2 2 -3.0
"""

        with tempfile.TemporaryDirectory() as tmpdir:
            file_path = self.create_test_file(tmpdir, content)

            reader = PolymatReader(file_path)

            polymers = list(reader.iterate_polymers())
            assert len(polymers) == 3

            # Check first polymer
            counts, fe, conc = polymers[0]
            assert np.array_equal(counts, np.array([1, 0]))
            assert fe == -1.0
            assert conc is None

            # Check second polymer
            counts, fe, conc = polymers[1]
            assert np.array_equal(counts, np.array([0, 1]))
            assert fe == -2.0
            assert conc is None

    def test_file_not_found(self):
        """Test error handling for non-existent file."""
        with pytest.raises(FileNotFoundError):
            PolymatReader("/nonexistent/file.tbnpolymat")

    def test_invalid_file_format(self):
        """Test handling of invalid file format."""
        content = """# No valid data lines, just comments
# And no Number of monomers specified"""

        with tempfile.TemporaryDirectory() as tmpdir:
            file_path = self.create_test_file(tmpdir, content)

            reader = PolymatReader(file_path)
            with pytest.raises(ValueError):
                reader.read()


class TestConvenienceFunctions:
    """Test convenience functions."""

    def test_load_polymat_file(self):
        """Test load_polymat_file convenience function."""
        content = """# TBN Polymer Matrix
# Number of polymers: 1
# Number of monomers: 2
# Columns: monomer_counts[1..2]
#
1 1
"""

        with tempfile.TemporaryDirectory() as tmpdir:
            file_path = os.path.join(tmpdir, "test.tbnpolymat")
            with open(file_path, "w") as f:
                f.write(content)

            data = load_polymat_file(file_path)
            assert data.n_polymers == 1
            assert data.n_monomers == 2

    def test_save_polymat_file(self):
        """Test save_polymat_file convenience function."""
        with tempfile.TemporaryDirectory() as tmpdir:
            file_path = os.path.join(tmpdir, "test.tbnpolymat")

            data = PolymatData(
                polymers=[np.array([1, 1])],
                n_monomers=2,
                n_polymers=1,
                has_free_energies=False,
                has_concentrations=False,
            )

            save_polymat_file(file_path, data)

            # Verify file was created
            assert os.path.exists(file_path)

            # Read it back
            loaded_data = load_polymat_file(file_path)
            assert loaded_data.n_polymers == 1
            assert loaded_data.n_monomers == 2

    def test_check_matrix_hash(self):
        """Test check_matrix_hash function."""
        content = """# TBN Polymer Matrix
# Number of polymers: 1
# Number of monomers: 2
\\MATRIX-HASH: hash123
# Columns: monomer_counts[1..2]
#
1 1
"""

        with tempfile.TemporaryDirectory() as tmpdir:
            file_path = os.path.join(tmpdir, "test.tbnpolymat")
            with open(file_path, "w") as f:
                f.write(content)

            # Test matching hash
            assert check_matrix_hash(file_path, "hash123")

            # Test non-matching hash
            assert not check_matrix_hash(file_path, "wrong_hash")

            # Test non-existent file
            assert not check_matrix_hash("/nonexistent/file.tbnpolymat", "hash123")


class TestRoundTrip:
    """Test that writing and reading preserve data."""

    def test_round_trip_all_fields(self):
        """Test complete round trip with all fields."""
        with tempfile.TemporaryDirectory() as tmpdir:
            file_path = os.path.join(tmpdir, "test.tbnpolymat")

            # Create test data
            original_data = PolymatData(
                polymers=[np.array([1, 0, 3]), np.array([2, 1, 0]), np.array([0, 0, 5])],
                n_monomers=3,
                n_polymers=3,
                matrix_hash="roundtrip_hash",
                free_energies=np.array([-1.5, -2.7, -4.2]),
                concentrations=np.array([1.2e-6, 3.4e-7, 5.6e-8]),
                concentration_units="uM",
                has_free_energies=True,
                has_concentrations=True,
            )

            # Write
            writer = PolymatWriter(file_path)
            writer.write(original_data)

            # Read
            reader = PolymatReader(file_path)
            loaded_data = reader.read()

            # Compare
            assert loaded_data.n_monomers == original_data.n_monomers
            assert loaded_data.n_polymers == original_data.n_polymers
            assert loaded_data.matrix_hash == original_data.matrix_hash
            assert loaded_data.concentration_units == original_data.concentration_units
            assert loaded_data.has_free_energies == original_data.has_free_energies
            assert loaded_data.has_concentrations == original_data.has_concentrations

            # Compare polymer data
            for i in range(original_data.n_polymers):
                assert np.array_equal(loaded_data.polymers[i], original_data.polymers[i])

            # Compare free energies
            assert np.allclose(loaded_data.free_energies, original_data.free_energies)

            # Compare concentrations
            assert np.allclose(loaded_data.concentrations, original_data.concentrations)
