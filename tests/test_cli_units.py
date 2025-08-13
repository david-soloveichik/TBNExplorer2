import unittest
import tempfile
import os
import argparse
from unittest.mock import patch, MagicMock
from tbnexplorer2.cli import main
from tbnexplorer2.units import VALID_UNITS


class TestCLIUnits(unittest.TestCase):
    """Test cases for CLI concentration units functionality."""
    
    def setUp(self):
        """Set up test fixtures."""
        # Create a temporary TBN file with concentrations
        self.temp_file = tempfile.NamedTemporaryFile(mode='w', suffix='.tbn', delete=False)
        self.temp_file.write("""# Test TBN file with concentrations
\\UNITS: nM
A: a a*, 100.0
B: b b*, 50.0
""")
        self.temp_file.close()
        self.temp_filename = self.temp_file.name
    
    def tearDown(self):
        """Clean up test fixtures."""
        if os.path.exists(self.temp_filename):
            os.unlink(self.temp_filename)
    
    @patch('tbnexplorer2.cli.TBN')
    @patch('tbnexplorer2.cli.TBNParser.parse_file')
    @patch('tbnexplorer2.cli.PolymerBasisComputer')
    @patch('tbnexplorer2.cli.NormalizRunner')
    def test_default_concentration_units(self, mock_normaliz, mock_computer, mock_parser, mock_tbn):
        """Test that default concentration units are nM."""
        # Mock the necessary components
        mock_parser.return_value = ([], {}, 'nM')
        mock_tbn_instance = MagicMock()
        mock_tbn_instance.check_star_limiting.return_value = (True, None)
        mock_tbn.return_value = mock_tbn_instance
        mock_normaliz_instance = MagicMock()
        mock_normaliz_instance.check_normaliz_available.return_value = True
        mock_normaliz.return_value = mock_normaliz_instance
        mock_computer_instance = MagicMock()
        mock_computer_instance.compute_polymer_basis.return_value = []
        mock_computer.return_value = mock_computer_instance
        
        # Test with no units argument (should default to nM)
        with patch('sys.argv', ['tbnexplorer2', self.temp_filename]):
            try:
                main()
            except SystemExit:
                pass  # Expected due to mocking
        
        # Verify TBN was created with default units
        mock_tbn.assert_called_once()
        args, kwargs = mock_tbn.call_args
        self.assertEqual(kwargs.get('concentration_units'), 'nM')
    
    @patch('tbnexplorer2.cli.TBN')
    @patch('tbnexplorer2.cli.TBNParser.parse_file')
    @patch('tbnexplorer2.cli.PolymerBasisComputer')
    @patch('tbnexplorer2.cli.NormalizRunner')
    def test_custom_concentration_units(self, mock_normaliz, mock_computer, mock_parser, mock_tbn):
        """Test setting custom concentration units."""
        # Mock the necessary components
        mock_parser.return_value = ([], {})
        mock_tbn_instance = MagicMock()
        mock_tbn_instance.check_star_limiting.return_value = (True, None)
        mock_tbn.return_value = mock_tbn_instance
        mock_normaliz_instance = MagicMock()
        mock_normaliz_instance.check_normaliz_available.return_value = True
        mock_normaliz.return_value = mock_normaliz_instance
        mock_computer_instance = MagicMock()
        mock_computer_instance.compute_polymer_basis.return_value = []
        mock_computer.return_value = mock_computer_instance
        
        # Test each valid unit
        for unit in VALID_UNITS:
            with self.subTest(unit=unit):
                # Mock parser to return this specific unit
                mock_parser.return_value = ([], {}, unit)
                
                with patch('sys.argv', ['tbnexplorer2', self.temp_filename]):
                    try:
                        main()
                    except SystemExit:
                        pass  # Expected due to mocking
                
                # Find the call with the specific unit
                tbn_calls = mock_tbn.call_args_list
                found_unit = False
                for call in tbn_calls:
                    args, kwargs = call
                    if kwargs.get('concentration_units') == unit:
                        found_unit = True
                        break
                
                self.assertTrue(found_unit, f"TBN was not called with unit {unit}")
                
                # Reset mock for next iteration
                mock_tbn.reset_mock()
    
    def test_file_without_units_no_concentrations_allowed(self):
        """Test that files without UNITS keyword cannot have concentrations."""
        # Create a temporary file without UNITS but with concentrations
        no_units_file = tempfile.NamedTemporaryFile(mode='w', suffix='.tbn', delete=False)
        no_units_file.write("""# Test TBN file without UNITS
A: a a*, 100.0
B: b b*, 50.0
""")
        no_units_file.close()
        
        try:
            # This should raise an error
            with patch('sys.argv', ['tbnexplorer2', no_units_file.name]):
                with self.assertRaises(SystemExit):
                    main()
        finally:
            os.unlink(no_units_file.name)


if __name__ == '__main__':
    unittest.main()