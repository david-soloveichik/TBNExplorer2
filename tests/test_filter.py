"""
Unit tests for the filter module.
"""

import pytest
import tempfile
import shutil
from pathlib import Path
import numpy as np
from tbnexplorer2.filter import PolymerFilter


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
        tbn_content = """B: b1 b2, 100
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
        assert filter.polymer_data['has_free_energies'] == True
        assert filter.polymer_data['has_concentrations'] == True
        assert filter.polymer_data['concentration_units'] == 'nM'
    
    def test_missing_polymat_file(self, temp_dir):
        """Test error handling when .tbnpolymat file is missing."""
        tbn_file = temp_dir / "missing.tbn"
        tbn_file.write_text("B: b1 b2, 100")
        
        with pytest.raises(FileNotFoundError, match="Cannot find polymer matrix file"):
            PolymerFilter(str(tbn_file))
    
    def test_load_polymat_data(self, sample_tbn_file, sample_polymat_file):
        """Test loading polymer data from .tbnpolymat file."""
        filter = PolymerFilter(str(sample_tbn_file))
        
        assert len(filter.polymer_data['polymers']) == 13
        assert filter.polymer_data['free_energies'] is not None
        assert len(filter.polymer_data['free_energies']) == 13
        assert filter.polymer_data['concentrations'] is not None
        assert len(filter.polymer_data['concentrations']) == 13
        
        # Check first polymer
        first_polymer = filter.polymer_data['polymers'][0]
        np.testing.assert_array_equal(first_polymer, [0, 0, 0, 0, 0, 0, 1, 1])
        assert filter.polymer_data['free_energies'][0] == -2.0
        assert np.isclose(filter.polymer_data['concentrations'][0], 99.9)
    
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
    
    def test_percent_limit_filtering(self, sample_tbn_file, sample_polymat_file):
        """Test filtering with percent limit."""
        filter = PolymerFilter(str(sample_tbn_file))
        
        # Get all polymers first
        all_results = filter.filter_by_monomers([])  # Empty list should match all
        
        # Actually, empty list won't match properly. Let's filter for everything with B
        all_with_B = filter.filter_by_monomers(['B'])
        
        # Now filter with 10% limit
        filtered_results = filter.filter_by_monomers(['B'], percent_limit=10.0)
        
        # Should have fewer results
        assert len(filtered_results) < len(all_with_B)
        
        # Check that all remaining have high concentration
        total_conc = np.sum(filter.polymer_data['concentrations'])
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
        assert "Free energy:" in output
    
    def test_format_output_with_percent_limit(self, sample_tbn_file, sample_polymat_file):
        """Test output formatting with percent limit."""
        filter = PolymerFilter(str(sample_tbn_file))
        
        results = filter.filter_by_monomers(['B'], percent_limit=5.0)
        output = filter.format_output(results, ['B'], percent_limit=5.0)
        
        assert "# Percent limit: 5.0%" in output
    
    def test_polymat_without_concentrations(self, temp_dir):
        """Test handling .tbnpolymat file without concentrations."""
        # Create TBN file
        tbn_file = temp_dir / "no_conc.tbn"
        tbn_file.write_text("a b, 0\na* b*, 0")
        
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
        
        assert filter.polymer_data['has_concentrations'] == False
        assert filter.polymer_data['concentrations'] is None
        
        # Should still be able to filter
        results = filter.filter_by_monomers(['a b'])
        assert len(results) == 1
        assert results[0][3] is None  # No concentration
    
    def test_polymat_without_free_energies(self, temp_dir):
        """Test handling .tbnpolymat file without free energies."""
        # Create TBN file
        tbn_file = temp_dir / "no_fe.tbn"
        tbn_file.write_text("a b, 0\na* b*, 0")
        
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
        
        assert filter.polymer_data['has_free_energies'] == False
        assert filter.polymer_data['free_energies'] is None
        assert filter.polymer_data['has_concentrations'] == False  # No conc without FE
        
        # Should still be able to filter
        results = filter.filter_by_monomers(['a b'])
        assert len(results) == 1
        assert results[0][2] is None  # No free energy
        assert results[0][3] is None  # No concentration