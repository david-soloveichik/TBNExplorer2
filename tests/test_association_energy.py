"""Tests for association energy penalty calculations."""

import math
import unittest

from tbnexplorer2.polymer_basis import (
    _bimolecular,
    _celcius_to_kelvin,
    _water_density_mol_per_L,
    compute_assoc_energy_penalty,
)


class TestAssociationEnergy(unittest.TestCase):
    """Test association energy penalty calculations."""

    def test_celcius_to_kelvin(self):
        """Test temperature conversion."""
        self.assertAlmostEqual(_celcius_to_kelvin(0), 273.15)
        self.assertAlmostEqual(_celcius_to_kelvin(37), 310.15)
        self.assertAlmostEqual(_celcius_to_kelvin(100), 373.15)
        self.assertAlmostEqual(_celcius_to_kelvin(-273.15), 0)

    def test_water_density(self):
        """Test water density calculation."""
        # At 37°C, water density should be around 55.3 mol/L
        density_37 = _water_density_mol_per_L(37.0)
        self.assertAlmostEqual(density_37, 55.3, delta=0.5)

        # At 25°C, water density should be around 55.5 mol/L
        density_25 = _water_density_mol_per_L(25.0)
        self.assertAlmostEqual(density_25, 55.5, delta=0.5)

        # At 4°C (maximum density), should be around 55.6 mol/L
        density_4 = _water_density_mol_per_L(4.0)
        self.assertAlmostEqual(density_4, 55.6, delta=0.5)

    def test_bimolecular(self):
        """Test bimolecular association term calculation."""
        # Test with zero G and H - should only have the water density term
        KB = 0.001987204259
        temp_c = 37.0
        temp_k = 310.15
        water_density = _water_density_mol_per_L(temp_c)

        expected = -KB * temp_k * math.log(water_density)
        result = _bimolecular(temp_c, 0.0, 0.0)
        self.assertAlmostEqual(result, expected, places=6)

        # Test with non-zero G and H
        G = 5.0
        H = 3.0
        expected = (G - H) * temp_k / 310.15 + H - KB * temp_k * math.log(water_density)
        result = _bimolecular(temp_c, G, H)
        self.assertAlmostEqual(result, expected, places=6)

    def test_assoc_energy_penalty_single_monomer(self):
        """Test that single monomer has zero association penalty."""
        penalty = compute_assoc_energy_penalty(1, 37.0, 5.0, 3.0)
        self.assertEqual(penalty, 0.0)

    def test_assoc_energy_penalty_dimer(self):
        """Test association penalty for a dimer."""
        temp_c = 37.0
        G = 5.0
        H = 3.0

        # For 2 monomers, penalty = bimolecular * (2 - 1) = bimolecular
        expected = _bimolecular(temp_c, G, H)
        result = compute_assoc_energy_penalty(2, temp_c, G, H)
        self.assertAlmostEqual(result, expected, places=6)

    def test_assoc_energy_penalty_polymer(self):
        """Test association penalty for larger polymers."""
        temp_c = 37.0
        G = 5.0
        H = 3.0
        bimol = _bimolecular(temp_c, G, H)

        # Test for various polymer sizes
        for n_monomers in [3, 5, 10]:
            expected = bimol * (n_monomers - 1)
            result = compute_assoc_energy_penalty(n_monomers, temp_c, G, H)
            self.assertAlmostEqual(result, expected, places=6)

    def test_assoc_energy_penalty_temperature_dependence(self):
        """Test that association penalty varies with temperature."""
        G = 5.0
        H = 3.0
        n_monomers = 5

        penalty_25 = compute_assoc_energy_penalty(n_monomers, 25.0, G, H)
        penalty_37 = compute_assoc_energy_penalty(n_monomers, 37.0, G, H)
        penalty_50 = compute_assoc_energy_penalty(n_monomers, 50.0, G, H)

        # Penalties should be different at different temperatures
        self.assertNotAlmostEqual(penalty_25, penalty_37, places=3)
        self.assertNotAlmostEqual(penalty_37, penalty_50, places=3)

    def test_assoc_energy_penalty_zero_parameters(self):
        """Test association penalty with zero G and H parameters."""
        temp_c = 37.0
        n_monomers = 5

        # With G=0 and H=0, should still have water density contribution
        penalty = compute_assoc_energy_penalty(n_monomers, temp_c, 0.0, 0.0)

        # Should be negative (favorable due to water density term)
        self.assertLess(penalty, 0)

        # Check exact value
        KB = 0.001987204259
        temp_k = 310.15
        water_density = _water_density_mol_per_L(temp_c)
        expected = -KB * temp_k * math.log(water_density) * (n_monomers - 1)
        self.assertAlmostEqual(penalty, expected, places=6)


if __name__ == "__main__":
    unittest.main()
