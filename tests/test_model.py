import pytest
import numpy as np
from tbnexplorer2.model import BindingSite, Monomer, TBN


class TestBindingSite:
    
    def test_binding_site_creation(self):
        """Test creating binding sites."""
        site1 = BindingSite("a", False)
        assert site1.name == "a"
        assert site1.is_star == False
        assert str(site1) == "a"
        
        site2 = BindingSite("b", True)
        assert site2.name == "b"
        assert site2.is_star == True
        assert str(site2) == "b*"
    
    def test_binding_site_equality(self):
        """Test binding site equality."""
        site1 = BindingSite("a", False)
        site2 = BindingSite("a", False)
        site3 = BindingSite("a", True)
        site4 = BindingSite("b", False)
        
        assert site1 == site2
        assert site1 != site3
        assert site1 != site4
    
    def test_binding_site_hash(self):
        """Test that binding sites can be used in sets."""
        sites = {
            BindingSite("a", False),
            BindingSite("a", True),
            BindingSite("b", False),
            BindingSite("a", False),  # Duplicate
        }
        assert len(sites) == 3


class TestMonomer:
    
    def test_monomer_creation(self):
        """Test creating monomers."""
        sites = [BindingSite("a", False), BindingSite("b", True)]
        monomer = Monomer("mol1", sites, 100.0, "mol1: a b*")
        
        assert monomer.name == "mol1"
        assert len(monomer.binding_sites) == 2
        assert monomer.concentration == 100.0
        assert monomer.original_line == "mol1: a b*"
    
    def test_monomer_to_vector(self):
        """Test converting monomer to vector representation."""
        sites = [
            BindingSite("a", False),
            BindingSite("a", True),
            BindingSite("b", False),
            BindingSite("b", True),
            BindingSite("c", False),
        ]
        monomer = Monomer(None, sites, None, "a a* b b* c")
        
        binding_site_index = {"a": 0, "b": 1, "c": 2}
        
        vector = monomer.to_vector(binding_site_index)
        
        # a: +1, a*: -1 => 0
        # b: +1, b*: -1 => 0
        # c: +1 => 1
        expected = np.array([0, 0, 1])
        np.testing.assert_array_equal(vector, expected)
    
    def test_monomer_to_vector_unbalanced(self):
        """Test monomer vector with unbalanced binding sites."""
        sites = [
            BindingSite("a", False),
            BindingSite("a", False),
            BindingSite("b", True),
        ]
        monomer = Monomer(None, sites, None, "a a b*")
        
        binding_site_index = {"a": 0, "b": 1}
        
        vector = monomer.to_vector(binding_site_index)
        
        # a: +2, b*: -1
        expected = np.array([2, -1])
        np.testing.assert_array_equal(vector, expected)
    
    def test_monomer_string_representation(self):
        """Test monomer string representation."""
        sites = [BindingSite("a", False), BindingSite("b", True)]
        
        # Named monomer
        monomer1 = Monomer("mol1", sites, None, "mol1: a b*")
        assert str(monomer1) == "mol1: a b*"
        
        # Unnamed monomer
        monomer2 = Monomer(None, sites, None, "a b*")
        assert str(monomer2) == "a b*"


class TestTBN:
    
    def create_simple_tbn(self):
        """Helper to create a simple TBN for testing."""
        sites1 = [BindingSite("a", False), BindingSite("b", True)]
        sites2 = [BindingSite("a", True), BindingSite("b", False)]
        
        monomers = [
            Monomer("M1", sites1, 100.0, "M1: a b*"),
            Monomer("M2", sites2, 50.0, "M2: a* b"),
        ]
        
        binding_site_index = {"a": 0, "b": 1}
        
        return TBN(monomers, binding_site_index, concentration_units='nM')
    
    def test_tbn_creation(self):
        """Test creating a TBN."""
        tbn = self.create_simple_tbn()
        
        assert len(tbn.monomers) == 2
        assert len(tbn.binding_site_index) == 2
    
    def test_matrix_A(self):
        """Test matrix A construction."""
        tbn = self.create_simple_tbn()
        
        A = tbn.matrix_A
        
        # Matrix A should be 2x2 (2 binding sites, 2 monomers)
        assert A.shape == (2, 2)
        
        # M1: a (+1), b* (-1)
        # M2: a* (-1), b (+1)
        expected = np.array([
            [1, -1],  # a row: M1 has a, M2 has a*
            [-1, 1],  # b row: M1 has b*, M2 has b
        ])
        np.testing.assert_array_equal(A, expected)
    
    def test_concentrations(self):
        """Test concentration vector in Molar units."""
        tbn = self.create_simple_tbn()
        
        # The concentrations property now returns values in Molar
        # With default nM units: 100 nM = 1e-7 M, 50 nM = 5e-8 M
        concentrations = tbn.concentrations
        expected = np.array([1e-7, 5e-8])
        np.testing.assert_array_almost_equal(concentrations, expected)
        
        # Test original units property
        original_concentrations = tbn.concentrations_original_units
        expected_original = np.array([100.0, 50.0])
        np.testing.assert_array_equal(original_concentrations, expected_original)
    
    def test_star_limiting_valid(self):
        """Test star-limiting check for valid TBN."""
        sites1 = [BindingSite("a", False), BindingSite("b", False)]
        sites2 = [BindingSite("a", True), BindingSite("b", True)]
        
        monomers = [
            Monomer(None, sites1, 100.0, "a b"),
            Monomer(None, sites2, 100.0, "a* b*"),
        ]
        
        binding_site_index = {"a": 0, "b": 1}
        tbn = TBN(monomers, binding_site_index)
        
        is_valid, error_msg = tbn.check_star_limiting()
        
        assert is_valid == True
        assert error_msg is None
    
    def test_star_limiting_invalid(self):
        """Test star-limiting check for invalid TBN."""
        sites1 = [BindingSite("a", True), BindingSite("a", True)]
        sites2 = [BindingSite("b", False)]
        
        monomers = [
            Monomer(None, sites1, 100.0, "a* a*"),
            Monomer(None, sites2, 50.0, "b"),
        ]
        
        binding_site_index = {"a": 0, "b": 1}
        tbn = TBN(monomers, binding_site_index)
        
        is_valid, error_msg = tbn.check_star_limiting()
        
        assert is_valid == False
        assert "not star-limited" in error_msg
        assert "a:" in error_msg
    
    def test_star_limiting_without_concentrations(self):
        """Test star-limiting with unit concentrations."""
        sites1 = [BindingSite("a", False), BindingSite("b", False)]
        sites2 = [BindingSite("a", True), BindingSite("b", True)]
        
        monomers = [
            Monomer(None, sites1, None, "a b"),
            Monomer(None, sites2, None, "a* b*"),
        ]
        
        binding_site_index = {"a": 0, "b": 1}
        tbn = TBN(monomers, binding_site_index)
        
        is_valid, error_msg = tbn.check_star_limiting()
        
        assert is_valid == True
        assert error_msg is None
    
    def test_augmented_matrix_no_singletons_needed(self):
        """Test augmented matrix when all binding sites have singletons."""
        # Create TBN with singleton star monomers for each binding site
        sites1 = [BindingSite("a", True)]
        sites2 = [BindingSite("b", True)]
        
        monomers = [
            Monomer(None, sites1, None, "a*"),
            Monomer(None, sites2, None, "b*"),
        ]
        
        binding_site_index = {"a": 0, "b": 1}
        tbn = TBN(monomers, binding_site_index)
        
        A_prime, n_original = tbn.get_augmented_matrix_for_polymer_basis()
        
        assert n_original == 2
        # No additional columns needed
        assert A_prime.shape == (2, 2)
    
    def test_augmented_matrix_with_singletons_needed(self):
        """Test augmented matrix when singleton monomers need to be added."""
        sites = [
            BindingSite("a", False),
            BindingSite("b", False),
            BindingSite("c", False),
        ]
        
        monomers = [
            Monomer(None, sites, None, "a b c"),
        ]
        
        binding_site_index = {"a": 0, "b": 1, "c": 2}
        tbn = TBN(monomers, binding_site_index)
        
        A_prime, n_original = tbn.get_augmented_matrix_for_polymer_basis()
        
        assert n_original == 1
        # Should add 3 singleton columns
        assert A_prime.shape == (3, 4)
        
        # Check that singleton columns were added correctly
        # Original column plus three singleton columns
        expected = np.array([
            [1, -1, 0, 0],   # a row
            [1, 0, -1, 0],   # b row
            [1, 0, 0, -1],   # c row
        ])
        np.testing.assert_array_equal(A_prime, expected)