"""
Tests for matrix hashing and polymer basis caching functionality.
"""

import pytest
import numpy as np
import tempfile
import os
from pathlib import Path

from tbnexplorer2.model import TBN, Monomer, BindingSite
from tbnexplorer2.polymer_basis import PolymerBasisComputer, Polymer
from tbnexplorer2.normaliz import NormalizRunner


class TestMatrixHashing:
    """Test matrix hash computation."""
    
    def test_compute_matrix_hash_consistency(self):
        """Test that hash computation is consistent."""
        # Create a simple TBN
        binding_sites = [
            BindingSite("a", False),
            BindingSite("a", True),
            BindingSite("b", False),
            BindingSite("b", True)
        ]
        
        monomers = [
            Monomer("M1", [binding_sites[0], binding_sites[1]], None, "a a*"),
            Monomer("M2", [binding_sites[2], binding_sites[3]], None, "b b*")
        ]
        
        binding_site_index = {"a": 0, "a*": 1, "b": 2, "b*": 3}
        
        tbn1 = TBN(monomers, binding_site_index)
        tbn2 = TBN(monomers, binding_site_index)
        
        # Hash should be the same for identical TBNs
        hash1 = tbn1.compute_matrix_hash()
        hash2 = tbn2.compute_matrix_hash()
        
        assert hash1 == hash2
        assert len(hash1) == 64  # SHA256 produces 64 hex characters
    
    def test_compute_matrix_hash_different_for_different_matrices(self):
        """Test that different matrices produce different hashes."""
        binding_site_index = {"a": 0, "a*": 1, "b": 2, "b*": 3}
        
        # First TBN
        monomers1 = [
            Monomer("M1", [BindingSite("a", False), BindingSite("a", True)], None, "a a*"),
            Monomer("M2", [BindingSite("b", False), BindingSite("b", True)], None, "b b*")
        ]
        tbn1 = TBN(monomers1, binding_site_index)
        
        # Second TBN with different monomers
        monomers2 = [
            Monomer("M1", [BindingSite("a", False), BindingSite("b", True)], None, "a b*"),
            Monomer("M2", [BindingSite("b", False), BindingSite("a", True)], None, "b a*")
        ]
        tbn2 = TBN(monomers2, binding_site_index)
        
        hash1 = tbn1.compute_matrix_hash()
        hash2 = tbn2.compute_matrix_hash()
        
        assert hash1 != hash2


class TestPolymerBasisCaching:
    """Test polymer basis caching functionality."""
    
    def create_sample_tbn(self):
        """Create a sample TBN for testing."""
        binding_sites = [
            BindingSite("a", False),
            BindingSite("a", True),
            BindingSite("b", False),
            BindingSite("b", True)
        ]
        
        monomers = [
            Monomer("M1", [binding_sites[0], binding_sites[1]], 100.0, "a a*, 100"),
            Monomer("M2", [binding_sites[2], binding_sites[3]], 50.0, "b b*, 50")
        ]
        
        binding_site_index = {"a": 0, "a*": 1, "b": 2, "b*": 3}
        return TBN(monomers, binding_site_index, "nM")
    
    def create_sample_polymat_file(self, temp_dir, matrix_hash, polymers_data):
        """Create a sample .tbnpolymat file."""
        polymat_file = temp_dir / "test.tbnpolymat"
        
        with open(polymat_file, 'w') as f:
            f.write("# TBN Polymer Matrix\n")
            f.write("# Number of polymers: 2\n")
            f.write("# Number of monomers: 2\n")
            f.write(f"\\MATRIX-HASH: {matrix_hash}\n")
            f.write("# Concentration units: nM\n")
            f.write("# Columns: monomer_counts[1..2] free_energy concentration\n")
            f.write("#\n")
            
            for polymer_data in polymers_data:
                f.write(" ".join(map(str, polymer_data)) + "\n")
        
        return polymat_file
    
    def test_load_cached_polymer_basis_hash_match(self):
        """Test loading cached polymer basis when hash matches."""
        tbn = self.create_sample_tbn()
        computer = PolymerBasisComputer(tbn, None)
        
        # Get the actual hash from the TBN
        matrix_hash = tbn.compute_matrix_hash()
        
        with tempfile.TemporaryDirectory() as temp_dir:
            temp_path = Path(temp_dir)
            
            # Create a sample .tbnpolymat file with matching hash
            polymers_data = [
                [1, 0, -1.0, 95.5],  # 1 M1, 0 M2, free_energy, concentration
                [0, 1, -1.0, 47.3]   # 0 M1, 1 M2, free_energy, concentration
            ]
            
            polymat_file = self.create_sample_polymat_file(temp_path, matrix_hash, polymers_data)
            
            # Load cached polymer basis
            cached_polymers = computer.load_cached_polymer_basis(str(polymat_file))
            
            assert cached_polymers is not None
            assert len(cached_polymers) == 2
            
            # Check first polymer
            assert np.array_equal(cached_polymers[0].monomer_counts, np.array([1, 0]))
            assert np.array_equal(cached_polymers[1].monomer_counts, np.array([0, 1]))
    
    def test_load_cached_polymer_basis_hash_mismatch(self):
        """Test that cached polymer basis is not loaded when hash doesn't match."""
        tbn = self.create_sample_tbn()
        computer = PolymerBasisComputer(tbn, None)
        
        with tempfile.TemporaryDirectory() as temp_dir:
            temp_path = Path(temp_dir)
            
            # Create a sample .tbnpolymat file with wrong hash
            wrong_hash = "abcd1234" * 8  # 64 character wrong hash
            polymers_data = [
                [1, 0, -1.0, 95.5],
                [0, 1, -1.0, 47.3]
            ]
            
            polymat_file = self.create_sample_polymat_file(temp_path, wrong_hash, polymers_data)
            
            # Load cached polymer basis
            cached_polymers = computer.load_cached_polymer_basis(str(polymat_file))
            
            assert cached_polymers is None
    
    def test_load_cached_polymer_basis_no_hash(self):
        """Test that cached polymer basis is not loaded when no hash is present."""
        tbn = self.create_sample_tbn()
        computer = PolymerBasisComputer(tbn, None)
        
        with tempfile.TemporaryDirectory() as temp_dir:
            temp_path = Path(temp_dir)
            polymat_file = temp_path / "test.tbnpolymat"
            
            # Create a .tbnpolymat file without hash
            with open(polymat_file, 'w') as f:
                f.write("# TBN Polymer Matrix\n")
                f.write("# Number of polymers: 1\n")
                f.write("# Number of monomers: 2\n")
                f.write("# Columns: monomer_counts[1..2]\n")
                f.write("#\n")
                f.write("1 0\n")
            
            # Load cached polymer basis
            cached_polymers = computer.load_cached_polymer_basis(str(polymat_file))
            
            assert cached_polymers is None
    
    def test_load_cached_polymer_basis_file_not_exists(self):
        """Test that None is returned when file doesn't exist."""
        tbn = self.create_sample_tbn()
        computer = PolymerBasisComputer(tbn, None)
        
        # Try to load from non-existent file
        cached_polymers = computer.load_cached_polymer_basis("/non/existent/file.tbnpolymat")
        
        assert cached_polymers is None
    
    def test_save_tbnpolymat_includes_hash(self):
        """Test that saving .tbnpolymat file includes the matrix hash."""
        tbn = self.create_sample_tbn()
        computer = PolymerBasisComputer(tbn, None)
        
        # Create sample polymers
        polymers = [
            Polymer(np.array([1, 0]), tbn.monomers, tbn),
            Polymer(np.array([0, 1]), tbn.monomers, tbn)
        ]
        
        with tempfile.TemporaryDirectory() as temp_dir:
            temp_path = Path(temp_dir)
            polymat_file = temp_path / "test.tbnpolymat"
            
            # Save the file
            computer.save_tbnpolymat(
                polymers,
                str(polymat_file),
                compute_free_energies=False,
                compute_concentrations=False
            )
            
            # Read the file and check for hash
            with open(polymat_file, 'r') as f:
                content = f.read()
            
            expected_hash = tbn.compute_matrix_hash()
            assert f"\\MATRIX-HASH: {expected_hash}" in content