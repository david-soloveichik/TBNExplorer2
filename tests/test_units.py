import unittest

import numpy as np

from tbnexplorer2.units import (
    UNIT_TO_MOLAR,
    VALID_UNITS,
    convert_concentration,
    from_molar,
    get_unit_display_name,
    to_molar,
    validate_unit,
)


class TestUnits(unittest.TestCase):
    """Test cases for concentration unit conversion utilities."""

    def test_valid_units_constants(self):
        """Test that the valid units constants are correctly defined."""
        expected_units = ["pM", "nM", "uM", "mM", "M"]
        self.assertEqual(set(VALID_UNITS), set(expected_units))

        # Test conversion factors
        self.assertEqual(UNIT_TO_MOLAR["pM"], 1e-12)
        self.assertEqual(UNIT_TO_MOLAR["nM"], 1e-9)
        self.assertEqual(UNIT_TO_MOLAR["uM"], 1e-6)
        self.assertEqual(UNIT_TO_MOLAR["mM"], 1e-3)
        self.assertEqual(UNIT_TO_MOLAR["M"], 1.0)

    def test_validate_unit(self):
        """Test unit validation."""
        # Valid units should pass
        for unit in VALID_UNITS:
            self.assertEqual(validate_unit(unit), unit)

        # Invalid units should raise ValueError
        invalid_units = ["pm", "nm", "um", "mm", "m", "L", "uL", "invalid"]
        for unit in invalid_units:
            with self.assertRaises(ValueError) as context:
                validate_unit(unit)
            self.assertIn("Invalid concentration unit", str(context.exception))
            self.assertIn("Supported units:", str(context.exception))

    def test_to_molar_single_values(self):
        """Test conversion of single values to Molar."""
        # Test each unit
        self.assertEqual(to_molar(1.0, "M"), 1.0)
        self.assertEqual(to_molar(1000.0, "mM"), 1.0)
        self.assertEqual(to_molar(1e6, "uM"), 1.0)
        self.assertEqual(to_molar(1e9, "nM"), 1.0)
        self.assertEqual(to_molar(1e12, "pM"), 1.0)

        # Test zero
        self.assertEqual(to_molar(0.0, "nM"), 0.0)

        # Test different values
        self.assertAlmostEqual(to_molar(100.0, "nM"), 1e-7)
        self.assertAlmostEqual(to_molar(50.0, "uM"), 5e-5)

    def test_to_molar_arrays(self):
        """Test conversion of arrays to Molar."""
        values = np.array([100.0, 50.0, 25.0])
        expected = np.array([1e-7, 5e-8, 2.5e-8])
        result = to_molar(values, "nM")
        np.testing.assert_array_almost_equal(result, expected)

    def test_from_molar_single_values(self):
        """Test conversion from Molar to other units."""
        # Test each unit
        self.assertEqual(from_molar(1.0, "M"), 1.0)
        self.assertAlmostEqual(from_molar(1.0, "mM"), 1000.0)
        self.assertAlmostEqual(from_molar(1.0, "uM"), 1e6)
        self.assertAlmostEqual(from_molar(1.0, "nM"), 1e9, places=5)
        self.assertAlmostEqual(from_molar(1.0, "pM"), 1e12, places=5)

        # Test zero
        self.assertEqual(from_molar(0.0, "nM"), 0.0)

        # Test different values
        self.assertAlmostEqual(from_molar(1e-7, "nM"), 100.0)
        self.assertAlmostEqual(from_molar(5e-5, "uM"), 50.0)

    def test_from_molar_arrays(self):
        """Test conversion from Molar arrays to other units."""
        values = np.array([1e-7, 5e-8, 2.5e-8])
        expected = np.array([100.0, 50.0, 25.0])
        result = from_molar(values, "nM")
        np.testing.assert_array_almost_equal(result, expected)

    def test_convert_concentration_same_unit(self):
        """Test conversion between same units (should be no-op)."""
        value = 100.0
        for unit in VALID_UNITS:
            result = convert_concentration(value, unit, unit)
            self.assertEqual(result, value)

    def test_convert_concentration_different_units(self):
        """Test conversion between different units."""
        # Test specific conversions
        self.assertAlmostEqual(convert_concentration(1.0, "M", "mM"), 1000.0)
        self.assertAlmostEqual(convert_concentration(1000.0, "mM", "M"), 1.0)
        self.assertAlmostEqual(convert_concentration(1.0, "uM", "nM"), 1000.0)
        self.assertAlmostEqual(convert_concentration(1000.0, "nM", "uM"), 1.0)
        self.assertAlmostEqual(convert_concentration(1.0, "nM", "pM"), 1000.0)

        # Test roundtrip conversion
        original = 100.0
        converted = convert_concentration(original, "nM", "uM")
        roundtrip = convert_concentration(converted, "uM", "nM")
        self.assertAlmostEqual(roundtrip, original)

    def test_convert_concentration_arrays(self):
        """Test conversion of concentration arrays."""
        values = np.array([1.0, 2.0, 3.0])
        expected = np.array([1000.0, 2000.0, 3000.0])
        result = convert_concentration(values, "uM", "nM")
        np.testing.assert_array_almost_equal(result, expected)

    def test_get_unit_display_name(self):
        """Test getting display names for units."""
        expected_names = {
            "pM": "picoMolar (pM)",
            "nM": "nanoMolar (nM)",
            "uM": "microMolar (uM)",
            "mM": "milliMolar (mM)",
            "M": "Molar (M)",
        }

        for unit, expected in expected_names.items():
            self.assertEqual(get_unit_display_name(unit), expected)

        # Test invalid unit
        with self.assertRaises(ValueError):
            get_unit_display_name("invalid")

    def test_edge_cases(self):
        """Test edge cases and special values."""
        # Very small values
        small_value = 1e-15
        result = to_molar(small_value, "pM")
        expected = small_value * 1e-12
        self.assertAlmostEqual(result, expected)

        # Very large values
        large_value = 1e15
        result = from_molar(large_value, "pM")
        expected = large_value * 1e12
        self.assertAlmostEqual(result, expected)

        # Negative values (should work mathematically)
        negative_value = -100.0
        result = to_molar(negative_value, "nM")
        self.assertAlmostEqual(result, -1e-7)


if __name__ == "__main__":
    unittest.main()
