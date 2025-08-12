#!/usr/bin/env python3
"""Tests for the CLI module."""

import unittest
import tempfile
import os
import sys
from pathlib import Path
from unittest.mock import patch, MagicMock
import argparse

# Add parent directory to path to import tbnexplorer2
sys.path.insert(0, str(Path(__file__).parent.parent))

from tbnexplorer2.cli import main


class TestCLIOutputFileLocation(unittest.TestCase):
    """Test that output files are created in the correct location."""
    
    def setUp(self):
        """Set up test fixtures."""
        # Create a temporary directory for tests
        self.test_dir = tempfile.mkdtemp()
        self.subdir = Path(self.test_dir) / "subdir"
        self.subdir.mkdir()
        
        # Create a simple valid TBN file content
        self.tbn_content = """# Test TBN file
monomer1: a b
monomer2: a* b*
"""
    
    def tearDown(self):
        """Clean up test fixtures."""
        import shutil
        shutil.rmtree(self.test_dir)
    
    def test_default_output_same_directory_as_input(self):
        """Test that default output file is created in the same directory as input."""
        # Create test TBN file in subdirectory
        input_file = self.subdir / "test.tbn"
        input_file.write_text(self.tbn_content)
        
        # Expected output file path
        expected_output = self.subdir / "test-polymer-basis.txt"
        
        # Mock the necessary components to avoid actual computation
        with patch('tbnexplorer2.cli.TBNParser') as mock_parser, \
             patch('tbnexplorer2.cli.TBN') as mock_tbn, \
             patch('tbnexplorer2.cli.NormalizRunner') as mock_normaliz, \
             patch('tbnexplorer2.cli.PolymerBasisComputer') as mock_computer:
            
            # Set up mocks
            mock_parser.parse_file.return_value = ([], {}, None)
            mock_tbn_instance = MagicMock()
            mock_tbn_instance.check_star_limiting.return_value = (True, None)
            mock_tbn_instance.matrix_A = MagicMock(shape=(0, 0))
            mock_tbn.return_value = mock_tbn_instance
            
            mock_normaliz_instance = MagicMock()
            mock_normaliz_instance.check_normaliz_available.return_value = True
            mock_normaliz.return_value = mock_normaliz_instance
            
            mock_computer_instance = MagicMock()
            mock_computer_instance.compute_polymer_basis.return_value = []
            mock_computer.return_value = mock_computer_instance
            
            # Run the CLI with the test file and user-friendly flag
            test_args = ['tbnexplorer2', str(input_file), '--user-friendly-polymer-basis']
            with patch('sys.argv', test_args):
                try:
                    main()
                except SystemExit:
                    pass  # main() calls sys.exit(0) on success
            
            # Verify that save_polymer_basis was called with the correct path
            mock_computer_instance.save_polymer_basis.assert_called_once()
            actual_output = mock_computer_instance.save_polymer_basis.call_args[0][1]
            self.assertEqual(actual_output, str(expected_output))
    
    def test_explicit_output_file_path_respected(self):
        """Test that explicitly specified output file path is used."""
        # Create test TBN file
        input_file = self.subdir / "test.tbn"
        input_file.write_text(self.tbn_content)
        
        # Explicit output file in different location
        explicit_output = Path(self.test_dir) / "custom-output.txt"
        
        # Mock the necessary components
        with patch('tbnexplorer2.cli.TBNParser') as mock_parser, \
             patch('tbnexplorer2.cli.TBN') as mock_tbn, \
             patch('tbnexplorer2.cli.NormalizRunner') as mock_normaliz, \
             patch('tbnexplorer2.cli.PolymerBasisComputer') as mock_computer:
            
            # Set up mocks
            mock_parser.parse_file.return_value = ([], {}, None)
            mock_tbn_instance = MagicMock()
            mock_tbn_instance.check_star_limiting.return_value = (True, None)
            mock_tbn_instance.matrix_A = MagicMock(shape=(0, 0))
            mock_tbn.return_value = mock_tbn_instance
            
            mock_normaliz_instance = MagicMock()
            mock_normaliz_instance.check_normaliz_available.return_value = True
            mock_normaliz.return_value = mock_normaliz_instance
            
            mock_computer_instance = MagicMock()
            mock_computer_instance.compute_polymer_basis.return_value = []
            mock_computer.return_value = mock_computer_instance
            
            # Run the CLI with explicit output file and user-friendly flag
            test_args = ['tbnexplorer2', str(input_file), '--user-friendly-polymer-basis', '--output', str(explicit_output)]
            with patch('sys.argv', test_args):
                try:
                    main()
                except SystemExit:
                    pass
            
            # Verify that save_polymer_basis was called with the explicit path
            mock_computer_instance.save_polymer_basis.assert_called_once()
            actual_output = mock_computer_instance.save_polymer_basis.call_args[0][1]
            self.assertEqual(actual_output, str(explicit_output))
    
    def test_current_directory_input_file(self):
        """Test that files in current directory work correctly."""
        # Create test TBN file in test directory
        input_file = Path(self.test_dir) / "current.tbn"
        input_file.write_text(self.tbn_content)
        
        # Expected output file path (same directory)
        expected_output = Path(self.test_dir) / "current-polymer-basis.txt"
        
        # Mock the necessary components
        with patch('tbnexplorer2.cli.TBNParser') as mock_parser, \
             patch('tbnexplorer2.cli.TBN') as mock_tbn, \
             patch('tbnexplorer2.cli.NormalizRunner') as mock_normaliz, \
             patch('tbnexplorer2.cli.PolymerBasisComputer') as mock_computer:
            
            # Set up mocks
            mock_parser.parse_file.return_value = ([], {}, None)
            mock_tbn_instance = MagicMock()
            mock_tbn_instance.check_star_limiting.return_value = (True, None)
            mock_tbn_instance.matrix_A = MagicMock(shape=(0, 0))
            mock_tbn.return_value = mock_tbn_instance
            
            mock_normaliz_instance = MagicMock()
            mock_normaliz_instance.check_normaliz_available.return_value = True
            mock_normaliz.return_value = mock_normaliz_instance
            
            mock_computer_instance = MagicMock()
            mock_computer_instance.compute_polymer_basis.return_value = []
            mock_computer.return_value = mock_computer_instance
            
            # Run the CLI with user-friendly flag
            test_args = ['tbnexplorer2', str(input_file), '--user-friendly-polymer-basis']
            with patch('sys.argv', test_args):
                try:
                    main()
                except SystemExit:
                    pass
            
            # Verify output path
            mock_computer_instance.save_polymer_basis.assert_called_once()
            actual_output = mock_computer_instance.save_polymer_basis.call_args[0][1]
            self.assertEqual(actual_output, str(expected_output))


class TestUserFriendlyPolymerBasisFlag(unittest.TestCase):
    """Test the --user-friendly-polymer-basis flag behavior."""
    
    def setUp(self):
        """Set up test fixtures."""
        # Create a temporary directory for tests
        self.test_dir = tempfile.mkdtemp()
        
        # Create a simple valid TBN file content
        self.tbn_content = """# Test TBN file
monomer1: a b
monomer2: a* b*
"""
        # Create test TBN file
        self.input_file = Path(self.test_dir) / "test.tbn"
        self.input_file.write_text(self.tbn_content)
    
    def tearDown(self):
        """Clean up test fixtures."""
        import shutil
        shutil.rmtree(self.test_dir)
    
    def test_user_friendly_flag_saves_basis_file(self):
        """Test that --user-friendly-polymer-basis flag causes polymer basis file to be saved."""
        with patch('tbnexplorer2.cli.TBNParser') as mock_parser, \
             patch('tbnexplorer2.cli.TBN') as mock_tbn, \
             patch('tbnexplorer2.cli.NormalizRunner') as mock_normaliz, \
             patch('tbnexplorer2.cli.PolymerBasisComputer') as mock_computer:
            
            # Set up mocks
            mock_parser.parse_file.return_value = ([], {}, None)
            mock_tbn_instance = MagicMock()
            mock_tbn_instance.check_star_limiting.return_value = (True, None)
            mock_tbn_instance.matrix_A = MagicMock(shape=(0, 0))
            mock_tbn_instance.concentrations = None
            mock_tbn.return_value = mock_tbn_instance
            
            mock_normaliz_instance = MagicMock()
            mock_normaliz_instance.check_normaliz_available.return_value = True
            mock_normaliz.return_value = mock_normaliz_instance
            
            mock_computer_instance = MagicMock()
            mock_computer_instance.compute_polymer_basis.return_value = []
            mock_computer.return_value = mock_computer_instance
            
            # Run CLI with --user-friendly-polymer-basis flag
            test_args = ['tbnexplorer2', str(self.input_file), '--user-friendly-polymer-basis']
            with patch('sys.argv', test_args):
                try:
                    main()
                except SystemExit:
                    pass
            
            # Verify that save_polymer_basis was called (user-friendly file should be saved)
            mock_computer_instance.save_polymer_basis.assert_called_once()
            
            # Verify that save_tbnpolymat was also called (always saved)
            mock_computer_instance.save_tbnpolymat.assert_called_once()
    
    def test_no_user_friendly_flag_skips_basis_file(self):
        """Test that without --user-friendly-polymer-basis flag, polymer basis file is not saved."""
        with patch('tbnexplorer2.cli.TBNParser') as mock_parser, \
             patch('tbnexplorer2.cli.TBN') as mock_tbn, \
             patch('tbnexplorer2.cli.NormalizRunner') as mock_normaliz, \
             patch('tbnexplorer2.cli.PolymerBasisComputer') as mock_computer:
            
            # Set up mocks
            mock_parser.parse_file.return_value = ([], {}, None)
            mock_tbn_instance = MagicMock()
            mock_tbn_instance.check_star_limiting.return_value = (True, None)
            mock_tbn_instance.matrix_A = MagicMock(shape=(0, 0))
            mock_tbn_instance.concentrations = None
            mock_tbn.return_value = mock_tbn_instance
            
            mock_normaliz_instance = MagicMock()
            mock_normaliz_instance.check_normaliz_available.return_value = True
            mock_normaliz.return_value = mock_normaliz_instance
            
            mock_computer_instance = MagicMock()
            mock_computer_instance.compute_polymer_basis.return_value = []
            mock_computer.return_value = mock_computer_instance
            
            # Run CLI without --user-friendly-polymer-basis flag
            test_args = ['tbnexplorer2', str(self.input_file)]
            with patch('sys.argv', test_args):
                try:
                    main()
                except SystemExit:
                    pass
            
            # Verify that save_polymer_basis was NOT called (no user-friendly file)
            mock_computer_instance.save_polymer_basis.assert_not_called()
            
            # Verify that save_tbnpolymat was called (always saved)
            mock_computer_instance.save_tbnpolymat.assert_called_once()


if __name__ == '__main__':
    unittest.main()