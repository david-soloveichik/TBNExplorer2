"""
Unit tests for the filter module.
"""

import pytest
import tempfile
import shutil
from pathlib import Path
import numpy as np
from tbnexplorer2.filter import PolymerFilter, format_concentration_nicely


class TestFormatConcentrationNicely:
    """Test the nice concentration formatting function."""
    
    def test_zero_value(self):
        """Test formatting of zero concentration."""
        assert format_concentration_nicely(0, 'nM') == '0.00 nM'
        assert format_concentration_nicely(0, 'uM') == '0.00 uM'
    
    def test_very_small_values(self):
        """Test formatting of very small values (<0.001)."""
        assert format_concentration_nicely(1e-5, 'nM') == '1.00e-05 nM'
        assert format_concentration_nicely(5.23e-10, 'nM') == '5.23e-10 nM'
    
    def test_small_values(self):
        """Test formatting of small values (0.001-1)."""
        assert format_concentration_nicely(0.001, 'nM') == '0.0010 nM'
        assert format_concentration_nicely(0.0056, 'nM') == '0.0056 nM'
        assert format_concentration_nicely(0.023, 'nM') == '0.023 nM'
        assert format_concentration_nicely(0.12, 'nM') == '0.12 nM'
        assert format_concentration_nicely(0.5, 'nM') == '0.50 nM'
    
    def test_medium_values(self):
        """Test formatting of medium values (1-1000)."""
        assert format_concentration_nicely(1.5, 'nM') == '1.50 nM'
        assert format_concentration_nicely(9.99, 'nM') == '9.99 nM'
        assert format_concentration_nicely(15.3, 'nM') == '15.3 nM'
        assert format_concentration_nicely(99.9, 'nM') == '99.9 nM'
        assert format_concentration_nicely(150.5, 'nM') == '150.5 nM'
        assert format_concentration_nicely(999.9, 'nM') == '999.9 nM'
    
    def test_large_values(self):
        """Test formatting of large values (1000-10000)."""
        assert format_concentration_nicely(1234, 'nM') == '1234 nM'
        assert format_concentration_nicely(5678.9, 'nM') == '5679 nM'
        assert format_concentration_nicely(9999, 'nM') == '9999 nM'
    
    def test_very_large_values(self):
        """Test formatting of very large values (>10000)."""
        assert format_concentration_nicely(15000, 'nM') == '1.50e+04 nM'
        assert format_concentration_nicely(1.23e8, 'nM') == '1.23e+08 nM'
    
    def test_negative_values(self):
        """Test formatting of negative values (shouldn't happen but handle gracefully)."""
        assert format_concentration_nicely(-10.5, 'nM') == '-10.5 nM'
        assert format_concentration_nicely(-0.001, 'nM') == '-0.0010 nM'
    
    def test_different_units(self):
        """Test formatting works with different units."""
        assert format_concentration_nicely(100, 'pM') == '100.0 pM'
        assert format_concentration_nicely(100, 'uM') == '100.0 uM'
        assert format_concentration_nicely(100, 'mM') == '100.0 mM'
        assert format_concentration_nicely(100, 'M') == '100.0 M'


class TestPolymerFilter:
    """Test the PolymerFilter class."""
    
    @pytest.fixture
    def temp_dir(self):
        """Create a temporary directory for test files."""
        temp_dir = tempfile.mkdtemp()
        yield Path(temp_dir)
        shutil.rmtree(temp_dir)
    
    @pytest.fixture
    def sample_tbn_file(self, temp_dir):
        """Create a sample .tbn file for testing."""
        tbn_content = """UNITS: nM
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
        return tbn_file
    
    @pytest.fixture
    def sample_polymat_file(self, temp_dir, sample_tbn_file):
        """Create a sample .tbnpolymat file for testing."""
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
        return polymat_file
    
    def test_init_and_file_inference(self, sample_tbn_file, sample_polymat_file):
        """Test initialization and .tbnpolymat file inference."""
        filter = PolymerFilter(str(sample_tbn_file))
        
        assert filter.tbn_file == sample_tbn_file
        assert filter.polymat_file == sample_polymat_file
        assert len(filter.monomers) == 8
        assert filter.polymer_data.has_free_energies == True
        assert filter.polymer_data.has_concentrations == True
        assert filter.units == 'nM'
    
    def test_missing_polymat_file(self, temp_dir):
        """Test error handling when .tbnpolymat file is missing."""
        tbn_file = temp_dir / "missing.tbn"
        tbn_file.write_text("UNITS: nM\nB: b1 b2, 100")
        
        with pytest.raises(FileNotFoundError, match="Cannot find polymer matrix file"):
            PolymerFilter(str(tbn_file))
    
    def test_load_polymat_data(self, sample_tbn_file, sample_polymat_file):
        """Test loading polymer data from .tbnpolymat file."""
        filter = PolymerFilter(str(sample_tbn_file))
        
        assert len(filter.polymer_data.polymers) == 13
        assert filter.polymer_data.free_energies is not None
        assert len(filter.polymer_data.free_energies) == 13
        assert filter.polymer_data.concentrations is not None
        assert len(filter.polymer_data.concentrations) == 13
        
        # Check first polymer
        first_polymer = filter.polymer_data.polymers[0]
        np.testing.assert_array_equal(first_polymer, [0, 0, 0, 0, 0, 0, 1, 1])
        assert filter.polymer_data.free_energies[0] == -2.0
        assert np.isclose(filter.polymer_data.concentrations[0], 99.9)
    
    def test_filter_by_single_monomer(self, sample_tbn_file, sample_polymat_file):
        """Test filtering by a single monomer name."""
        filter = PolymerFilter(str(sample_tbn_file))
        
        # Filter for polymers containing "B"
        results = filter.filter_by_monomers(['B'])
        
        # Should find polymers with monomer index 0 (B)
        assert len(results) > 0
        
        # Check that all results contain monomer B
        for idx, counts, fe, conc in results:
            assert counts[0] > 0  # B is at index 0
    
    def test_filter_by_multiple_monomers(self, sample_tbn_file, sample_polymat_file):
        """Test filtering by multiple monomer names."""
        filter = PolymerFilter(str(sample_tbn_file))
        
        # Filter for polymers containing both "B" and "C"
        results = filter.filter_by_monomers(['B', 'C'])
        
        # Should find polymers with both monomers
        assert len(results) > 0
        
        for idx, counts, fe, conc in results:
            assert counts[0] > 0  # B is at index 0
            assert counts[7] > 0  # C is at index 7
    
    def test_filter_with_multiplicity(self, sample_tbn_file, sample_polymat_file):
        """Test filtering with monomer multiplicity (duplicates)."""
        filter = PolymerFilter(str(sample_tbn_file))
        
        # Filter for polymers containing at least 2 monomers with binding sites c1* c2*
        # This corresponds to monomer at index 6
        results = filter.filter_by_monomers(['c1* c2*', 'c1* c2*'])
        
        # Should only find polymers with count >= 2 for monomer 6
        for idx, counts, fe, conc in results:
            assert counts[6] >= 2
    
    def test_filter_nonexistent_monomer(self, sample_tbn_file, sample_polymat_file):
        """Test filtering with a non-existent monomer name."""
        filter = PolymerFilter(str(sample_tbn_file))
        
        # Filter for non-existent monomer
        results = filter.filter_by_monomers(['NonExistent'])
        
        # Should return empty list
        assert len(results) == 0
    
    def test_filter_empty_monomer_list(self, sample_tbn_file, sample_polymat_file):
        """Test filtering with empty monomer list (should return all polymers)."""
        filter = PolymerFilter(str(sample_tbn_file))
        
        # Filter with empty list
        results = filter.filter_by_monomers([])
        
        # Should return all 13 polymers
        assert len(results) == 13
        
        # All polymers should be included
        polymer_indices = set(idx for idx, _, _, _ in results)
        assert polymer_indices == set(range(13))
    
    def test_max_count_limit(self, sample_tbn_file, sample_polymat_file):
        """Test filtering with max_count limit."""
        filter = PolymerFilter(str(sample_tbn_file))
        
        # Get all polymers first
        all_results = filter.filter_by_monomers([])
        
        # Filter with max_count=5
        limited_results = filter.filter_by_monomers([], max_count=5)
        
        # Should have exactly 5 results
        assert len(limited_results) == 5
        
        # Should be the first 5 from the full list (after sorting)
        for i in range(5):
            assert limited_results[i][0] == all_results[i][0]  # Same polymer index
    
    def test_max_count_with_monomer_filter(self, sample_tbn_file, sample_polymat_file):
        """Test max_count combined with monomer filtering."""
        filter = PolymerFilter(str(sample_tbn_file))
        
        # Filter for B with max_count=2
        results = filter.filter_by_monomers(['B'], max_count=2)
        
        # Should have at most 2 results
        assert len(results) <= 2
        
        # All should contain B
        for idx, counts, fe, conc in results:
            assert counts[0] > 0  # B is at index 0
    
    def test_percent_limit_filtering(self, sample_tbn_file, sample_polymat_file):
        """Test filtering with percent limit."""
        filter = PolymerFilter(str(sample_tbn_file))
        
        # Get all polymers with B first
        all_with_B = filter.filter_by_monomers(['B'])
        
        # Now filter with 10% limit
        filtered_results = filter.filter_by_monomers(['B'], percent_limit=10.0)
        
        # Should have fewer results
        assert len(filtered_results) <= len(all_with_B)
        
        # Check that all remaining have high concentration
        total_conc = np.sum(filter.polymer_data.concentrations)
        for idx, counts, fe, conc in filtered_results:
            assert (conc / total_conc * 100) >= 10.0
    
    def test_sorting_by_concentration(self, sample_tbn_file, sample_polymat_file):
        """Test that results are sorted by decreasing concentration."""
        filter = PolymerFilter(str(sample_tbn_file))
        
        # Get polymers containing B
        results = filter.filter_by_monomers(['B'])
        
        # Check sorting
        concentrations = [conc for _, _, _, conc in results]
        assert concentrations == sorted(concentrations, reverse=True)
    
    def test_format_output(self, sample_tbn_file, sample_polymat_file):
        """Test output formatting."""
        filter = PolymerFilter(str(sample_tbn_file))
        
        # Filter for polymers containing B and C
        results = filter.filter_by_monomers(['B', 'C'])
        
        # Format output
        output = filter.format_output(results, ['B', 'C'])
        
        # Check output contains expected elements
        assert "# Filtered polymers containing: B C" in output
        assert "# Number of matching polymers:" in output
        assert "# Total concentration fraction:" in output
        assert "# Concentration units: nM" in output
        assert "# Polymer 1" in output
        assert "Concentration:" in output
        # Free energy is no longer included in output to save space
        assert "Free energy:" not in output
    
    def test_format_output_with_percent_limit(self, sample_tbn_file, sample_polymat_file):
        """Test output formatting with percent limit."""
        filter = PolymerFilter(str(sample_tbn_file))
        
        results = filter.filter_by_monomers(['B'], percent_limit=5.0)
        output = filter.format_output(results, ['B'], percent_limit=5.0)
        
        assert "# Percent limit: 5.0%" in output
    
    def test_format_output_with_max_count(self, sample_tbn_file, sample_polymat_file):
        """Test output formatting with max_count limit."""
        filter = PolymerFilter(str(sample_tbn_file))
        
        results = filter.filter_by_monomers(['B'], max_count=3)
        output = filter.format_output(results, ['B'], max_count=3)
        
        assert "# Maximum count limit: 3" in output
    
    def test_format_output_empty_monomers(self, sample_tbn_file, sample_polymat_file):
        """Test output formatting with empty monomer list."""
        filter = PolymerFilter(str(sample_tbn_file))
        
        results = filter.filter_by_monomers([])
        output = filter.format_output(results, [])
        
        assert "# All polymers" in output
        assert "# Filtered polymers containing:" not in output
    
    def test_nice_concentration_formatting_in_output(self, sample_tbn_file, sample_polymat_file):
        """Test that concentrations are nicely formatted in output."""
        filter = PolymerFilter(str(sample_tbn_file))
        
        # Get polymers with concentrations that span different ranges
        results = filter.filter_by_monomers(['C'])
        output = filter.format_output(results, ['C'])
        
        # Check for nice formatting of concentrations
        # From the sample data: 9.99e+01 should become "99.9 nM"
        assert "99.9 nM" in output
        # 1.22e-01 should become "0.12 nM"
        assert "0.12 nM" in output
        # Very small values should still use scientific notation
        assert "1.48e-04 nM" in output
    
    def test_polymat_without_concentrations(self, temp_dir):
        """Test handling .tbnpolymat file without concentrations."""
        # Create TBN file
        tbn_file = temp_dir / "no_conc.tbn"
        tbn_file.write_text("UNITS: nM\na b, 0\na* b*, 0")
        
        # Create polymat file without concentrations
        polymat_content = """# TBN Polymer Matrix
# Number of polymers: 2
# Number of monomers: 2
# Partial computation: no monomer concentrations provided
# Columns: monomer_counts[1..2] free_energy
#
1 0 -0.0
0 1 -0.0"""
        
        polymat_file = temp_dir / "no_conc.tbnpolymat"
        polymat_file.write_text(polymat_content)
        
        filter = PolymerFilter(str(tbn_file))
        
        assert filter.polymer_data.has_concentrations == False
        assert filter.polymer_data.concentrations is None
        
        # Should still be able to filter
        results = filter.filter_by_monomers(['a b'])
        assert len(results) == 1
        assert results[0][3] is None  # No concentration
    
    def test_polymat_without_free_energies(self, temp_dir):
        """Test handling .tbnpolymat file without free energies."""
        # Create TBN file
        tbn_file = temp_dir / "no_fe.tbn"
        tbn_file.write_text("UNITS: nM\na b, 0\na* b*, 0")
        
        # Create polymat file without free energies
        polymat_content = """# TBN Polymer Matrix
# Number of polymers: 2
# Number of monomers: 2
# Partial computation: free energies not computed
# Columns: monomer_counts[1..2]
#
1 0
0 1"""
        
        polymat_file = temp_dir / "no_fe.tbnpolymat"
        polymat_file.write_text(polymat_content)
        
        filter = PolymerFilter(str(tbn_file))
        
        assert filter.polymer_data.has_free_energies == False
        assert filter.polymer_data.free_energies is None
        assert filter.polymer_data.has_concentrations == False  # No conc without FE
        
        # Should still be able to filter
        results = filter.filter_by_monomers(['a b'])
        assert len(results) == 1
        assert results[0][2] is None  # No free energy
        assert results[0][3] is None  # No concentration