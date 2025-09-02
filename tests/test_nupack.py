"""Tests for NUPACK concentrations integration."""

import tempfile
from unittest.mock import Mock, patch

import numpy as np
import pytest

from tbnexplorer2.config import NUPACK_CONCENTRATIONS_PATH
from tbnexplorer2.model import TBN, BindingSite, Monomer
from tbnexplorer2.nupack import NupackRunner
from tbnexplorer2.polymer_basis import Polymer


@pytest.fixture
def sample_polymers():
    """Create sample polymers for testing."""
    # Create simple polymers with monomer counts using dummy monomers
    monomers = [Mock(spec=Monomer) for _ in range(3)]
    polymer1 = Polymer(np.array([2, 0, 1]), monomers)  # 2 of monomer 0, 1 of monomer 2
    polymer2 = Polymer(np.array([1, 1, 0]), monomers)  # 1 of monomer 0, 1 of monomer 1
    polymer3 = Polymer(np.array([0, 2, 1]), monomers)  # 2 of monomer 1, 1 of monomer 2
    return [polymer1, polymer2, polymer3]


@pytest.fixture
def sample_tbn():
    """Create a sample TBN with monomer concentrations."""
    monomers = [
        Monomer(
            name="A",
            binding_sites=[BindingSite("x", False), BindingSite("y", False)],
            concentration=1e-6,
            original_line="A: x y, 1e-6",
        ),  # 1 μM
        Monomer(
            name="B",
            binding_sites=[BindingSite("y", False), BindingSite("z", False)],
            concentration=2e-6,
            original_line="B: y z, 2e-6",
        ),  # 2 μM
        Monomer(
            name="C",
            binding_sites=[BindingSite("x", False), BindingSite("z", False)],
            concentration=0.5e-6,
            original_line="C: x z, 0.5e-6",
        ),  # 0.5 μM
    ]
    binding_site_index = {"x": 0, "y": 1, "z": 2}
    return TBN(monomers, binding_site_index, concentration_units="M")


class TestNupackRunner:
    """Test NupackRunner class."""

    def test_init(self):
        """Test NupackRunner initialization."""
        runner = NupackRunner("/path/to/nupack", temperature=25.0)
        assert runner.nupack_path == "/path/to/nupack"
        assert runner.temperature == 25.0

    def test_init_defaults(self):
        """Test NupackRunner with default values."""
        runner = NupackRunner()
        assert runner.nupack_path == NUPACK_CONCENTRATIONS_PATH  # Default from config
        assert runner.temperature == 37.0

    @patch("os.path.isfile")
    @patch("os.access")
    def test_check_nupack_available(self, mock_access, mock_isfile):
        """Test checking if NUPACK is available."""
        runner = NupackRunner("/path/to/nupack")

        # Test when NUPACK is available
        mock_isfile.return_value = True
        mock_access.return_value = True
        assert runner.check_nupack_available() is True

        # Test when file doesn't exist
        mock_isfile.return_value = False
        assert runner.check_nupack_available() is False

        # Test when file exists but not executable
        mock_isfile.return_value = True
        mock_access.return_value = False
        assert runner.check_nupack_available() is False

    def test_write_ocx_file(self, sample_polymers, tmp_path):
        """Test writing OCX file for NUPACK."""
        runner = NupackRunner()
        ocx_path = tmp_path / "test.ocx"

        # Mock compute_free_energy to return predictable values
        for i, polymer in enumerate(sample_polymers):
            polymer.compute_free_energy = Mock(return_value=-1.5 * (i + 1))

        runner._write_ocx_file(sample_polymers, str(ocx_path), deltaG=[0.0, 0.0], temperature=37.0)

        # Read and verify file contents
        with open(ocx_path) as f:
            lines = f.readlines()

        assert len(lines) == 3
        # Check format: id, 1, monomer counts, free energy
        assert lines[0].strip() == "1\t1\t2\t0\t1\t-1.5"
        assert lines[1].strip() == "2\t1\t1\t1\t0\t-3.0"
        assert lines[2].strip() == "3\t1\t0\t2\t1\t-4.5"

    def test_write_con_file(self, sample_tbn, tmp_path):
        """Test writing CON file for NUPACK."""
        runner = NupackRunner()
        con_path = tmp_path / "test.con"

        runner._write_con_file(sample_tbn, str(con_path))

        # Read and verify file contents
        with open(con_path) as f:
            lines = f.readlines()

        assert len(lines) == 3
        assert float(lines[0].strip()) == pytest.approx(1e-6)
        assert float(lines[1].strip()) == pytest.approx(2e-6)
        assert float(lines[2].strip()) == pytest.approx(0.5e-6)

    def test_parse_nupack_output(self, tmp_path):
        """Test parsing NUPACK .eq output file."""
        runner = NupackRunner()
        eq_path = tmp_path / "test.eq"

        # Create a sample .eq file
        eq_content = """1\t1\t2\t0\t1\t-1.5\t1.23e-7
2\t1\t1\t1\t0\t-3.0\t4.56e-8
3\t1\t0\t2\t1\t-4.5\t7.89e-9"""
        with open(eq_path, "w") as f:
            f.write(eq_content)

        concentrations = runner._parse_nupack_output(str(eq_path))

        assert len(concentrations) == 3
        assert concentrations[0] == pytest.approx(1.23e-7)
        assert concentrations[1] == pytest.approx(4.56e-8)
        assert concentrations[2] == pytest.approx(7.89e-9)

    def test_parse_nupack_output_invalid_format(self, tmp_path):
        """Test parsing invalid NUPACK output."""
        runner = NupackRunner()
        eq_path = tmp_path / "test.eq"

        # Create invalid .eq file
        with open(eq_path, "w") as f:
            f.write("invalid\n")

        with pytest.raises(RuntimeError, match="Invalid .eq file format"):
            runner._parse_nupack_output(str(eq_path))

    @patch("subprocess.run")
    def test_compute_equilibrium_concentrations(self, mock_run, sample_polymers, sample_tbn, tmp_path):
        """Test full equilibrium concentration computation."""
        runner = NupackRunner("/path/to/nupack", temperature=25.0)

        # Mock NUPACK availability check
        runner.check_nupack_available = Mock(return_value=True)

        # Mock subprocess to simulate successful NUPACK run
        mock_result = Mock()
        mock_result.returncode = 0
        mock_result.stderr = ""
        mock_run.return_value = mock_result

        # Mock compute_free_energy for polymers
        for i, polymer in enumerate(sample_polymers):
            polymer.compute_free_energy = Mock(return_value=-1.5 * (i + 1))

        # Create temporary .eq file that NUPACK would generate
        with tempfile.TemporaryDirectory() as temp_dir:
            # Mock the parse method to return concentrations
            runner._parse_nupack_output = Mock(return_value=np.array([1e-7, 2e-8, 3e-9]))

            concentrations = runner.compute_equilibrium_concentrations(
                sample_polymers, sample_tbn, output_dir=temp_dir, deltaG=[0.0, 0.0], temperature=25.0
            )

            # Verify subprocess was called with correct arguments
            assert mock_run.called
            call_args = mock_run.call_args[0][0]
            assert call_args[0] == "/path/to/nupack"
            assert call_args[1] == "-sort"
            assert call_args[2] == "0"
            assert call_args[3] == "-T"
            assert call_args[4] == "25.0"
            assert "nupack_input" in call_args[5]

            # Verify concentrations returned
            assert len(concentrations) == 3
            assert concentrations[0] == pytest.approx(1e-7)
            assert concentrations[1] == pytest.approx(2e-8)
            assert concentrations[2] == pytest.approx(3e-9)

    def test_compute_equilibrium_concentrations_no_concentrations(self, sample_polymers):
        """Test that computation fails without monomer concentrations."""
        runner = NupackRunner()
        # Create TBN without concentrations
        monomers = [
            Monomer("A", [BindingSite("x", False)], None, "A: x"),
            Monomer("B", [BindingSite("y", False)], None, "B: y"),
        ]
        binding_site_index = {"x": 0, "y": 1}
        tbn = TBN(monomers, binding_site_index)

        with pytest.raises(ValueError, match="Cannot compute equilibrium concentrations"):
            runner.compute_equilibrium_concentrations(sample_polymers, tbn)

    @patch("subprocess.run")
    def test_compute_equilibrium_concentrations_nupack_failure(self, mock_run, sample_polymers, sample_tbn):
        """Test handling of NUPACK execution failure."""
        runner = NupackRunner("/path/to/nupack")
        runner.check_nupack_available = Mock(return_value=True)

        # Mock subprocess to simulate NUPACK failure
        mock_result = Mock()
        mock_result.returncode = 1
        mock_result.stderr = "NUPACK error: invalid input"
        mock_run.return_value = mock_result

        # Mock compute_free_energy for polymers
        for polymer in sample_polymers:
            polymer.compute_free_energy = Mock(return_value=-1.0)

        with pytest.raises(RuntimeError, match="NUPACK failed: NUPACK error: invalid input"):
            runner.compute_equilibrium_concentrations(sample_polymers, sample_tbn)

    def test_compute_equilibrium_concentrations_nupack_not_found(self, sample_polymers, sample_tbn):
        """Test error when NUPACK is not available."""
        runner = NupackRunner("/nonexistent/nupack")
        runner.check_nupack_available = Mock(return_value=False)

        with pytest.raises(RuntimeError, match="NUPACK not found"):
            runner.compute_equilibrium_concentrations(sample_polymers, sample_tbn)
