"""Tests for constraints file functionality in tbnexplorer2-filter."""

import os
import subprocess
import sys
import tempfile
from pathlib import Path


class TestConstraintsFile:
    """Test constraints file functionality."""

    def setup_method(self):
        """Create temporary test files."""
        # Create a simple TBN file with units and concentrations
        with tempfile.NamedTemporaryFile(mode="w", suffix=".tbn", delete=False) as test_tbn:
            test_tbn.write("\\UNITS: nM\n")
            test_tbn.write("a b >M1, 100\n")
            test_tbn.write("a* b* >M2, 100\n")
            test_tbn.write("c d >M3, 100\n")
            test_tbn.write("c* d* >M4, 100\n")
            test_tbn.write("a a b b, 50\n")  # Unnamed monomer
            test_tbn.write("c c d d, 50\n")  # Another unnamed monomer
            self.test_tbn_name = test_tbn.name

        # Create corresponding polymat file
        with tempfile.NamedTemporaryFile(mode="w", suffix=".tbnpolymat", delete=False):
            base_name = Path(self.test_tbn_name).stem
            self.test_polymat_path = Path(self.test_tbn_name).parent / f"{base_name}.tbnpolymat"

        # Write polymat content
        with open(self.test_polymat_path, "w") as f:
            f.write("# TBN Polymer Matrix\n")
            f.write("# Number of polymers: 6\n")
            f.write("# Number of monomers: 6\n")
            f.write("\\MATRIX-HASH: test_hash\n")
            f.write("# Concentration units: nanoMolar (nM)\n")
            f.write("# Columns: monomer_counts[1..6] free_energy concentration\n")
            f.write("#\n")
            # Polymer 1: M1 + M2 dimer
            f.write("1 1 0 0 0 0 -2.0 50.0\n")
            # Polymer 2: M3 + M4 dimer
            f.write("0 0 1 1 0 0 -2.0 40.0\n")
            # Polymer 3: Just M1
            f.write("1 0 0 0 0 0 0.0 30.0\n")
            # Polymer 4: Just M3
            f.write("0 0 1 0 0 0 0.0 20.0\n")
            # Polymer 5: M1 + M3
            f.write("1 0 1 0 0 0 0.0 10.0\n")
            # Polymer 6: Unnamed monomer (a a b b)
            f.write("0 0 0 0 1 0 0.0 5.0\n")

    def teardown_method(self):
        """Clean up temporary files."""
        try:
            os.unlink(self.test_tbn_name)
            os.unlink(self.test_polymat_path)
        except Exception:
            pass

    def run_filter(self, args):
        """Run tbnexplorer2-filter with given arguments."""
        cmd = [sys.executable, "-m", "tbnexplorer2.filter_cli", *args]
        result = subprocess.run(cmd, capture_output=True, text=True)
        return result

    def test_contains_single_constraint(self):
        """Test CONTAINS constraint with single monomer."""
        # Create constraints file
        with tempfile.NamedTemporaryFile(mode="w", suffix=".txt", delete=False) as f:
            f.write("CONTAINS M1\n")
            constraints_file = f.name

        try:
            result = self.run_filter([self.test_tbn_name, "--constraints-file", constraints_file])
            assert result.returncode == 0
            assert "M1" in result.stdout
            # Should match polymers 1, 3, and 5
            assert "# Number of matching polymers: 3" in result.stdout
        finally:
            os.unlink(constraints_file)

    def test_contains_multiple_monomers(self):
        """Test CONTAINS constraint with multiple monomers."""
        with tempfile.NamedTemporaryFile(mode="w", suffix=".txt", delete=False) as f:
            f.write("CONTAINS M1 M2\n")
            constraints_file = f.name

        try:
            result = self.run_filter([self.test_tbn_name, "--constraints-file", constraints_file])
            assert result.returncode == 0
            # Should only match polymer 1 (M1 + M2 dimer)
            assert "# Number of matching polymers: 1" in result.stdout
        finally:
            os.unlink(constraints_file)

    def test_exactly_constraint(self):
        """Test EXACTLY constraint."""
        with tempfile.NamedTemporaryFile(mode="w", suffix=".txt", delete=False) as f:
            f.write("EXACTLY M1\n")
            constraints_file = f.name

        try:
            result = self.run_filter([self.test_tbn_name, "--constraints-file", constraints_file])
            assert result.returncode == 0
            # Should only match polymer 3 (just M1)
            assert "# Number of matching polymers: 1" in result.stdout
            assert "M1" in result.stdout
        finally:
            os.unlink(constraints_file)

    def test_or_logic_multiple_constraints(self):
        """Test OR logic with multiple constraints."""
        with tempfile.NamedTemporaryFile(mode="w", suffix=".txt", delete=False) as f:
            f.write("EXACTLY M1\n")
            f.write("EXACTLY M3\n")
            constraints_file = f.name

        try:
            result = self.run_filter([self.test_tbn_name, "--constraints-file", constraints_file])
            assert result.returncode == 0
            # Should match polymers 3 (just M1) and 4 (just M3)
            assert "# Number of matching polymers: 2" in result.stdout
        finally:
            os.unlink(constraints_file)

    def test_constraints_with_nonexistent_monomer(self):
        """Test constraints with monomer name that doesn't exist."""
        with tempfile.NamedTemporaryFile(mode="w", suffix=".txt", delete=False) as f:
            f.write("CONTAINS NonExistentMonomer\n")
            constraints_file = f.name

        try:
            result = self.run_filter([self.test_tbn_name, "--constraints-file", constraints_file])
            assert result.returncode == 0
            # Should match no polymers
            assert "# Number of matching polymers: 0" in result.stdout
        finally:
            os.unlink(constraints_file)

    def test_empty_constraints_file(self):
        """Test empty constraints file returns all polymers."""
        with tempfile.NamedTemporaryFile(mode="w", suffix=".txt", delete=False) as f:
            f.write("# Just comments\n")
            f.write("\n")
            constraints_file = f.name

        try:
            result = self.run_filter([self.test_tbn_name, "--constraints-file", constraints_file])
            assert result.returncode == 0
            # Should return all polymers
            assert "# Number of matching polymers: 6" in result.stdout
        finally:
            os.unlink(constraints_file)

    def test_invalid_constraint_type(self):
        """Test error handling for invalid constraint type."""
        with tempfile.NamedTemporaryFile(mode="w", suffix=".txt", delete=False) as f:
            f.write("INCLUDES M1\n")  # Invalid - should be CONTAINS
            constraints_file = f.name

        try:
            result = self.run_filter([self.test_tbn_name, "--constraints-file", constraints_file])
            assert result.returncode != 0
            assert "Invalid constraint type 'INCLUDES'" in result.stderr
        finally:
            os.unlink(constraints_file)

    def test_constraints_file_with_command_line_monomers(self):
        """Test error when specifying both constraints file and command line monomers."""
        with tempfile.NamedTemporaryFile(mode="w", suffix=".txt", delete=False) as f:
            f.write("CONTAINS M1\n")
            constraints_file = f.name

        try:
            result = self.run_filter([self.test_tbn_name, "M1", "--constraints-file", constraints_file])
            assert result.returncode != 0
            assert "Cannot specify monomer names on command line when using --constraints-file" in result.stderr
        finally:
            os.unlink(constraints_file)

    def test_nonexistent_constraints_file(self):
        """Test error handling for nonexistent constraints file."""
        result = self.run_filter([self.test_tbn_name, "--constraints-file", "nonexistent.txt"])
        assert result.returncode != 0
        assert "Constraints file 'nonexistent.txt' not found" in result.stderr

    def test_constraints_with_percent_limit(self):
        """Test constraints file with percent limit."""
        with tempfile.NamedTemporaryFile(mode="w", suffix=".txt", delete=False) as f:
            f.write("CONTAINS M1\n")
            constraints_file = f.name

        try:
            result = self.run_filter(
                [self.test_tbn_name, "--constraints-file", constraints_file, "--percent-limit", "20"]
            )
            assert result.returncode == 0
            # With percent limit, should only show polymers > 20% of total
            # Total concentration is 155, so > 31 nM
            # Only polymer 1 (50.0) should match
            assert "# Number of matching polymers: 1" in result.stdout
        finally:
            os.unlink(constraints_file)

    def test_constraints_with_num_limit(self):
        """Test constraints file with num limit."""
        with tempfile.NamedTemporaryFile(mode="w", suffix=".txt", delete=False) as f:
            f.write("CONTAINS M1\n")
            f.write("CONTAINS M3\n")
            constraints_file = f.name

        try:
            result = self.run_filter([self.test_tbn_name, "--constraints-file", constraints_file, "--num", "2"])
            assert result.returncode == 0
            # Should limit to 2 polymers
            assert "# Maximum count limit: 2" in result.stdout
            assert "# Number of matching polymers: 2" in result.stdout
        finally:
            os.unlink(constraints_file)
