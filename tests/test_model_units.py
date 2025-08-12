import unittest
import numpy as np
from tbnexplorer2.model import TBN, Monomer, BindingSite
from tbnexplorer2.units import to_molar


class TestTBNModelUnits(unittest.TestCase):
    """Test cases for TBN model concentration unit handling."""
    
    def setUp(self):
        """Set up test fixtures."""
        # Create test monomers with concentrations
        self.monomer1 = Monomer(
            name="A",
            binding_sites=[BindingSite("a", False), BindingSite("a", True)],
            concentration=100.0,
            original_line="A: a a*, 100.0"
        )
        
        self.monomer2 = Monomer(
            name="B", 
            binding_sites=[BindingSite("b", False), BindingSite("b", True)],
            concentration=50.0,
            original_line="B: b b*, 50.0"
        )
        
        self.monomers = [self.monomer1, self.monomer2]
        self.binding_site_index = {"a": 0, "b": 1}
    
    def test_no_concentration_units(self):
        """Test that concentration units default to None when not specified."""
        tbn = TBN(self.monomers, self.binding_site_index)
        self.assertIsNone(tbn.concentration_units)
    
    def test_custom_concentration_units(self):
        """Test setting custom concentration units."""
        for unit in ['pM', 'nM', 'uM', 'mM', 'M']:
            tbn = TBN(self.monomers, self.binding_site_index, concentration_units=unit)
            self.assertEqual(tbn.concentration_units, unit)
    
    def test_concentrations_original_units(self):
        """Test that original concentrations are preserved."""
        tbn = TBN(self.monomers, self.binding_site_index, concentration_units='nM')
        original_conc = tbn.concentrations_original_units
        
        self.assertIsNotNone(original_conc)
        expected = np.array([100.0, 50.0])
        np.testing.assert_array_equal(original_conc, expected)
    
    def test_concentrations_converted_to_molar(self):
        """Test that internal concentrations are converted to Molar."""
        # Test with nM units
        tbn = TBN(self.monomers, self.binding_site_index, concentration_units='nM')
        molar_conc = tbn.concentrations
        
        self.assertIsNotNone(molar_conc)
        # 100 nM = 100 * 1e-9 M = 1e-7 M
        # 50 nM = 50 * 1e-9 M = 5e-8 M
        expected = np.array([1e-7, 5e-8])
        np.testing.assert_array_almost_equal(molar_conc, expected)
    
    def test_concentrations_different_units(self):
        """Test concentration conversion with different input units."""
        test_cases = [
            ('pM', np.array([1e-10, 5e-11])),  # 100 pM, 50 pM -> Molar
            ('uM', np.array([1e-4, 5e-5])),   # 100 uM, 50 uM -> Molar
            ('mM', np.array([1e-1, 5e-2])),   # 100 mM, 50 mM -> Molar  
            ('M', np.array([100.0, 50.0])),   # 100 M, 50 M -> Molar
        ]
        
        for unit, expected_molar in test_cases:
            with self.subTest(unit=unit):
                tbn = TBN(self.monomers, self.binding_site_index, concentration_units=unit)
                molar_conc = tbn.concentrations
                np.testing.assert_array_almost_equal(molar_conc, expected_molar)
    
    def test_no_concentrations_returns_none(self):
        """Test that models without concentrations return None."""
        # Create monomers without concentrations
        monomer_no_conc = Monomer(
            name="A",
            binding_sites=[BindingSite("a", False)],
            concentration=None,
            original_line="A: a"
        )
        
        tbn = TBN([monomer_no_conc], {"a": 0})
        self.assertIsNone(tbn.concentrations)
        self.assertIsNone(tbn.concentrations_original_units)
    
    def test_mixed_concentrations_returns_none(self):
        """Test that models with mixed concentration specifications return None."""
        # One monomer with concentration, one without
        monomer_with_conc = Monomer(
            name="A",
            binding_sites=[BindingSite("a", False)],
            concentration=100.0,
            original_line="A: a, 100.0"
        )
        
        monomer_no_conc = Monomer(
            name="B", 
            binding_sites=[BindingSite("b", False)],
            concentration=None,
            original_line="B: b"
        )
        
        tbn = TBN([monomer_with_conc, monomer_no_conc], {"a": 0, "b": 1})
        self.assertIsNone(tbn.concentrations)
        self.assertIsNone(tbn.concentrations_original_units)
    
    def test_star_limiting_with_units(self):
        """Test that star limiting check works correctly with different units."""
        # Create a TBN that should be star-limited regardless of units
        # Both monomers have equal amounts of starred and unstarred sites
        tbn = TBN(self.monomers, self.binding_site_index, concentration_units='nM')
        
        is_valid, error_msg = tbn.check_star_limiting()
        self.assertTrue(is_valid)
        self.assertIsNone(error_msg)
        
        # Test with different units - should still be valid
        for unit in ['pM', 'uM', 'mM', 'M']:
            tbn_unit = TBN(self.monomers, self.binding_site_index, concentration_units=unit)
            is_valid, error_msg = tbn_unit.check_star_limiting()
            self.assertTrue(is_valid, f"Failed for unit {unit}")
            self.assertIsNone(error_msg)
    
    def test_concentration_caching(self):
        """Test that concentration calculations are cached properly."""
        tbn = TBN(self.monomers, self.binding_site_index, concentration_units='nM')
        
        # First call should compute and cache
        conc1 = tbn.concentrations
        original1 = tbn.concentrations_original_units
        
        # Second call should return cached values (same object)
        conc2 = tbn.concentrations
        original2 = tbn.concentrations_original_units
        
        self.assertIs(conc1, conc2)
        self.assertIs(original1, original2)


if __name__ == '__main__':
    unittest.main()