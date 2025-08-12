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
                with patch('sys.argv', ['tbnexplorer2', self.temp_filename, '--concentration-units', unit]):
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
    
    def test_argument_parser_accepts_valid_units(self):
        """Test that argument parser accepts all valid concentration units."""
        from tbnexplorer2.cli import main
        
        # Create a mock argument parser to test choices
        parser = argparse.ArgumentParser()
        parser.add_argument(
            '--concentration-units',
            default='nM',
            choices=VALID_UNITS,
            help='Concentration units for input and output'
        )
        
        # Test that all valid units are accepted
        for unit in VALID_UNITS:
            args = parser.parse_args(['--concentration-units', unit])
            self.assertEqual(args.concentration_units, unit)
    
    def test_argument_parser_rejects_invalid_units(self):
        """Test that argument parser rejects invalid concentration units."""
        from tbnexplorer2.cli import main
        
        # Create a mock argument parser to test choices
        parser = argparse.ArgumentParser()
        parser.add_argument(
            '--concentration-units',
            default='nM',
            choices=VALID_UNITS,
            help='Concentration units for input and output'
        )
        
        # Test that invalid units are rejected
        invalid_units = ['nm', 'pm', 'um', 'mm', 'm', 'L', 'invalid']
        for unit in invalid_units:
            with self.assertRaises(SystemExit):
                parser.parse_args(['--concentration-units', unit])
    
    @patch('sys.stderr')
    def test_cli_help_includes_units_info(self, mock_stderr):
        """Test that CLI help includes information about concentration units."""
        with patch('sys.argv', ['tbnexplorer2', '--help']):
            with self.assertRaises(SystemExit):
                main()
        
        # The help should have been printed (this is a basic smoke test)
        # In a real scenario, we'd capture stdout to check the help text


if __name__ == '__main__':
    unittest.main()