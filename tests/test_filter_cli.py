"""
Integration tests for the filter CLI.
"""

import shutil
import subprocess
import tempfile
from pathlib import Path

import pytest


class TestFilterCLI:
    """Test the tbnexplorer2-filter command line interface."""

    @pytest.fixture
    def temp_dir(self):
        """Create a temporary directory for test files."""
        temp_dir = tempfile.mkdtemp()
        yield Path(temp_dir)
        shutil.rmtree(temp_dir)

    @pytest.fixture
    def sample_files(self, temp_dir):
        """Create sample .tbn and .tbnpolymat files for testing."""
        # Create .tbn file
        tbn_content = """\\UNITS: nM
B: b1 b2, 100
a1* a2* b1* b2*, 100
a1 a2 b1 b2 c1, 100
a2* b1* b2* c1*, 100
a2 b1, 100
b2 c1 c2, 100
c1* c2*, 100
C: c1 c2, 100"""

        tbn_file = temp_dir / "test.tbn"
        tbn_file.write_text(tbn_content)

        # Create .tbnpolymat file
        polymat_content = """# TBN Polymer Matrix
# Number of polymers: 13
# Number of monomers: 8
# Concentration units: nM
# Columns: monomer_counts[1..8] free_energy concentration
#
0 0 0 0 0 0 1 1 -2.0 9.99e+01
0 0 0 1 1 1 0 0 -4.0 9.99e+01
0 1 1 0 0 0 0 0 -4.0 9.99e+01
1 0 0 0 0 0 0 0 -0.0 9.99e+01
0 0 0 0 0 0 0 1 -0.0 1.22e-01
0 0 0 0 0 1 1 0 -2.0 1.22e-01
1 1 1 1 1 0 0 0 -8.0 1.22e-01
0 0 0 0 0 1 0 0 -0.0 1.48e-04
0 1 1 1 1 1 1 0 -10.0 1.48e-04
1 0 0 1 1 0 0 1 -4.0 1.48e-04
0 0 1 1 0 0 0 0 -4.0 2.35e-11
0 0 0 0 1 0 0 0 -0.0 1.74e-11
0 0 1 0 0 0 0 0 -0.0 1.10e-38"""

        polymat_file = temp_dir / "test.tbnpolymat"
        polymat_file.write_text(polymat_content)

        return tbn_file, polymat_file

    def test_filter_single_monomer(self, sample_files):
        """Test filtering by a single monomer."""
        tbn_file, _ = sample_files

        result = subprocess.run(
            ["python", "-m", "tbnexplorer2.filter_cli", str(tbn_file), "B"], capture_output=True, text=True
        )

        assert result.returncode == 0
        assert "# Filtered polymers containing: B" in result.stdout
        assert "# Number of matching polymers:" in result.stdout
        assert "# Total concentration fraction:" in result.stdout
        assert "# Polymer 1" in result.stdout

    def test_filter_multiple_monomers(self, sample_files):
        """Test filtering by multiple monomers."""
        tbn_file, _ = sample_files

        result = subprocess.run(
            ["python", "-m", "tbnexplorer2.filter_cli", str(tbn_file), "B", "C"], capture_output=True, text=True
        )

        assert result.returncode == 0
        assert "# Filtered polymers containing: B C" in result.stdout
        assert "# Polymer" in result.stdout

    def test_filter_with_multiplicity(self, sample_files):
        """Test filtering with monomer multiplicity."""
        tbn_file, _ = sample_files

        result = subprocess.run(
            ["python", "-m", "tbnexplorer2.filter_cli", str(tbn_file), "B", "B"], capture_output=True, text=True
        )

        assert result.returncode == 0
        assert "# Filtered polymers containing: B B" in result.stdout

    def test_filter_with_percent_limit(self, sample_files):
        """Test filtering with percent limit."""
        tbn_file, _ = sample_files

        result = subprocess.run(
            ["python", "-m", "tbnexplorer2.filter_cli", str(tbn_file), "B", "--percent-limit", "10"],
            capture_output=True,
            text=True,
        )

        assert result.returncode == 0
        assert "# Percent limit: 10.0%" in result.stdout
        assert "# Filtered polymers containing: B" in result.stdout

    def test_filter_with_percent_limit_short_form(self, sample_files):
        """Test filtering with percent limit using short form -p."""
        tbn_file, _ = sample_files

        result = subprocess.run(
            ["python", "-m", "tbnexplorer2.filter_cli", str(tbn_file), "B", "-p", "10"], capture_output=True, text=True
        )

        assert result.returncode == 0
        assert "# Percent limit: 10.0%" in result.stdout
        assert "# Filtered polymers containing: B" in result.stdout

    def test_filter_with_num_limit(self, sample_files):
        """Test filtering with --num limit."""
        tbn_file, _ = sample_files

        result = subprocess.run(
            ["python", "-m", "tbnexplorer2.filter_cli", str(tbn_file), "B", "--num", "3"],
            capture_output=True,
            text=True,
        )

        assert result.returncode == 0
        assert "# Maximum count limit: 3" in result.stdout
        assert "# Filtered polymers containing: B" in result.stdout

    def test_filter_with_num_limit_short_form(self, sample_files):
        """Test filtering with -n limit."""
        tbn_file, _ = sample_files

        result = subprocess.run(
            ["python", "-m", "tbnexplorer2.filter_cli", str(tbn_file), "B", "-n", "2"], capture_output=True, text=True
        )

        assert result.returncode == 0
        assert "# Maximum count limit: 2" in result.stdout
        assert "# Filtered polymers containing: B" in result.stdout

    def test_filter_no_monomers(self, sample_files):
        """Test filtering with no monomers specified (should return all polymers)."""
        tbn_file, _ = sample_files

        result = subprocess.run(
            ["python", "-m", "tbnexplorer2.filter_cli", str(tbn_file)], capture_output=True, text=True
        )

        assert result.returncode == 0
        assert "# All polymers" in result.stdout
        assert "# Number of matching polymers: 13" in result.stdout  # Should match total from sample data

    def test_filter_no_monomers_with_num_limit(self, sample_files):
        """Test filtering with no monomers but with num limit."""
        tbn_file, _ = sample_files

        result = subprocess.run(
            ["python", "-m", "tbnexplorer2.filter_cli", str(tbn_file), "--num", "5"], capture_output=True, text=True
        )

        assert result.returncode == 0
        assert "# All polymers" in result.stdout
        assert "# Maximum count limit: 5" in result.stdout
        assert "# Number of matching polymers: 5" in result.stdout

    def test_invalid_percent_limit(self, sample_files):
        """Test error handling for invalid percent limit."""
        tbn_file, _ = sample_files

        # Test negative percent
        result = subprocess.run(
            ["python", "-m", "tbnexplorer2.filter_cli", str(tbn_file), "B", "--percent-limit", "-5"],
            capture_output=True,
            text=True,
        )

        assert result.returncode != 0
        assert "Error: --percent-limit must be between 0 and 100" in result.stderr

        # Test percent > 100
        result = subprocess.run(
            ["python", "-m", "tbnexplorer2.filter_cli", str(tbn_file), "B", "--percent-limit", "150"],
            capture_output=True,
            text=True,
        )

        assert result.returncode != 0
        assert "Error: --percent-limit must be between 0 and 100" in result.stderr

    def test_invalid_num_limit(self, sample_files):
        """Test error handling for invalid num limit."""
        tbn_file, _ = sample_files

        # Test negative num
        result = subprocess.run(
            ["python", "-m", "tbnexplorer2.filter_cli", str(tbn_file), "B", "--num", "0"],
            capture_output=True,
            text=True,
        )

        assert result.returncode != 0
        assert "Error: --num must be at least 1" in result.stderr

    def test_missing_tbn_file(self):
        """Test error handling for missing .tbn file."""
        result = subprocess.run(
            ["python", "-m", "tbnexplorer2.filter_cli", "nonexistent.tbn", "B"], capture_output=True, text=True
        )

        assert result.returncode != 0
        assert "Error: Input file 'nonexistent.tbn' not found" in result.stderr

    def test_missing_polymat_file(self, temp_dir):
        """Test error handling for missing .tbnpolymat file."""
        # Create only .tbn file without corresponding .tbnpolymat
        tbn_file = temp_dir / "incomplete.tbn"
        tbn_file.write_text("B: b1 b2, 100")

        result = subprocess.run(
            ["python", "-m", "tbnexplorer2.filter_cli", str(tbn_file), "B"], capture_output=True, text=True
        )

        assert result.returncode != 0
        assert "Cannot find polymer matrix file" in result.stderr

    def test_tbn_file_without_units(self, temp_dir):
        """Test error handling for .tbn file without UNITS keyword."""
        # Create .tbn file without UNITS keyword
        tbn_content = """# TBN file without UNITS
B: b1 b2
a1* a2* b1* b2*
a1 a2 b1 b2 c1"""

        tbn_file = temp_dir / "no_units.tbn"
        tbn_file.write_text(tbn_content)

        # Also create a dummy .tbnpolymat file to avoid the missing file error
        polymat_content = """# TBN Polymer Matrix
# Number of polymers: 1
# Number of monomers: 3
# Columns: monomer_counts[1..3]
#
1 0 0"""

        polymat_file = temp_dir / "no_units.tbnpolymat"
        polymat_file.write_text(polymat_content)

        result = subprocess.run(
            ["python", "-m", "tbnexplorer2.filter_cli", str(tbn_file), "B"], capture_output=True, text=True
        )

        assert result.returncode != 0
        assert "tbnexplorer2-filter requires a .tbn file with UNITS keyword and concentrations" in result.stderr
        assert "does not have UNITS specified" in result.stderr

    def test_nonexistent_monomer(self, sample_files):
        """Test filtering with non-existent monomer name."""
        tbn_file, _ = sample_files

        result = subprocess.run(
            ["python", "-m", "tbnexplorer2.filter_cli", str(tbn_file), "NonExistent"], capture_output=True, text=True
        )

        assert result.returncode == 0
        assert "# Number of matching polymers: 0" in result.stdout

    def test_help_message(self):
        """Test that help message is displayed correctly."""
        result = subprocess.run(["python", "-m", "tbnexplorer2.filter_cli", "--help"], capture_output=True, text=True)

        assert result.returncode == 0
        assert "Filter polymers from .tbnpolymat files by monomer names" in result.stdout
        assert "--percent-limit" in result.stdout
        assert "--num" in result.stdout
        assert "Examples:" in result.stdout

    def test_real_example_file(self):
        """Test with real example files if they exist."""
        # Check if the my_inputs directory exists
        my_inputs_dir = Path("my_inputs")
        if not my_inputs_dir.exists():
            pytest.skip("my_inputs directory not found")

        tbn_file = my_inputs_dir / "and_gate.tbn"
        polymat_file = my_inputs_dir / "and_gate.tbnpolymat"

        if not tbn_file.exists() or not polymat_file.exists():
            pytest.skip("Example files not found")

        # Test filtering
        result = subprocess.run(
            ["python", "-m", "tbnexplorer2.filter_cli", str(tbn_file), "B", "C"], capture_output=True, text=True
        )

        assert result.returncode == 0
        assert "# Filtered polymers containing: B C" in result.stdout
        assert "# Concentration:" in result.stdout
        # Free energy is no longer included in output to save space
        assert "Free energy:" not in result.stdout
