import os
import tempfile
from unittest.mock import MagicMock, Mock, patch

import numpy as np
import numpy.testing  # Import early to avoid subprocess issues during tests
import pytest

from tbnexplorer2.coffee import COFFEERunner
from tbnexplorer2.model import TBN
from tbnexplorer2.polymer_basis import Polymer


class TestCOFFEERunner:
    def test_init_custom_path(self):
        """Test COFFEERunner initialization with custom path."""
        custom_path = "/custom/path/to/coffee"
        runner = COFFEERunner(coffee_path=custom_path)
        assert runner.coffee_path == custom_path
        assert runner.temperature == 37.0  # Default temperature

    def test_init_with_temperature(self):
        """Test COFFEERunner initialization with temperature."""
        runner = COFFEERunner(temperature=25.0)
        assert runner.temperature == 25.0

    def test_init_default_temperature(self):
        """Test COFFEERunner uses default temperature of 37."""
        runner = COFFEERunner()
        assert runner.temperature == 37.0

    def test_check_coffee_available_exists(self):
        """Test check_coffee_available when file exists and is executable."""
        runner = COFFEERunner("/path/to/coffee")
        with patch("os.path.isfile", return_value=True), patch("os.access", return_value=True):
            assert runner.check_coffee_available() is True

    def test_check_coffee_available_not_exists(self):
        """Test check_coffee_available when file doesn't exist."""
        runner = COFFEERunner("/nonexistent/coffee")
        with patch("os.path.isfile", return_value=False):
            assert runner.check_coffee_available() is False

    def test_check_coffee_available_not_executable(self):
        """Test check_coffee_available when file exists but not executable."""
        runner = COFFEERunner("/path/to/coffee")
        with patch("os.path.isfile", return_value=True), patch("os.access", return_value=False):
            assert runner.check_coffee_available() is False

    def test_compute_equilibrium_no_concentrations(self):
        """Test compute_equilibrium_concentrations raises error without concentrations."""
        runner = COFFEERunner()
        tbn = MagicMock(spec=TBN)
        tbn.concentrations = None
        polymers = []

        with pytest.raises(ValueError, match="Cannot compute equilibrium concentrations"):
            runner.compute_equilibrium_concentrations(polymers, tbn)

    def test_compute_equilibrium_coffee_not_available(self):
        """Test compute_equilibrium_concentrations when COFFEE not available."""
        runner = COFFEERunner("/nonexistent/coffee")
        tbn = MagicMock(spec=TBN)
        tbn.concentrations = np.array([1.0, 2.0])
        polymers = []

        with patch.object(runner, "check_coffee_available", return_value=False):
            with pytest.raises(RuntimeError, match="COFFEE not found"):
                runner.compute_equilibrium_concentrations(polymers, tbn)

    def test_write_cfe_file(self):
        """Test _write_cfe_file method."""
        runner = COFFEERunner()

        # Create mock polymers
        polymer1 = Mock(spec=Polymer)
        polymer1.monomer_counts = np.array([1, 0, 1])
        polymer1.compute_free_energy = Mock(return_value=-2.5)

        polymer2 = Mock(spec=Polymer)
        polymer2.monomer_counts = np.array([0, 2, 1])
        polymer2.compute_free_energy = Mock(return_value=-3.0)

        polymers = [polymer1, polymer2]

        with tempfile.NamedTemporaryFile(mode="w", delete=False, suffix=".cfe") as f:
            filename = f.name

        try:
            runner._write_cfe_file(polymers, filename)

            # Read and verify the file
            with open(filename) as f:
                lines = f.readlines()

            assert len(lines) == 2
            assert lines[0].strip() == "1 0 1 -2.5"
            assert lines[1].strip() == "0 2 1 -3.0"
        finally:
            os.unlink(filename)

    def test_write_con_file(self):
        """Test _write_con_file method."""
        runner = COFFEERunner()
        tbn = Mock(spec=TBN)

        # tbn.concentrations should already return values in Molar units
        # These correspond to 100 nM, 50 nM, 25 nM after conversion to Molar
        tbn.concentrations = np.array([1e-7, 5e-8, 2.5e-8])  # Already in Molar
        tbn.concentration_units = "nM"

        with tempfile.NamedTemporaryFile(mode="w", delete=False, suffix=".con") as f:
            filename = f.name

        try:
            runner._write_con_file(tbn, filename)

            # Read and verify the file
            with open(filename) as f:
                lines = f.readlines()

            assert len(lines) == 3
            # Values should be written as-is since they're already in Molar
            assert float(lines[0].strip()) == pytest.approx(1e-7)
            assert float(lines[1].strip()) == pytest.approx(5e-8)
            assert float(lines[2].strip()) == pytest.approx(2.5e-8)
        finally:
            os.unlink(filename)

    def test_parse_coffee_output(self):
        """Test _parse_coffee_output method."""
        runner = COFFEERunner()

        # Create sample output file
        output_content = """1.23e-7
4.56e-8
0.00e0
7.89e-9
"""
        with tempfile.NamedTemporaryFile(mode="w", delete=False, suffix=".txt") as f:
            f.write(output_content)
            filename = f.name

        try:
            concentrations = runner._parse_coffee_output(filename)

            assert len(concentrations) == 4
            assert concentrations[0] == pytest.approx(1.23e-7)
            assert concentrations[1] == pytest.approx(4.56e-8)
            assert concentrations[2] == pytest.approx(0.0)
            assert concentrations[3] == pytest.approx(7.89e-9)
        finally:
            os.unlink(filename)

    def test_compute_equilibrium_concentrations_success(self):
        """Test successful equilibrium concentration computation."""
        runner = COFFEERunner("/path/to/coffee")

        # Mock TBN
        tbn = Mock(spec=TBN)
        tbn.concentrations = np.array([1e-7, 5e-8])  # Already in Molar (equivalent to 100 nM, 50 nM)
        tbn.concentration_units = "nM"

        # Mock polymers
        polymer1 = Mock(spec=Polymer)
        polymer1.monomer_counts = np.array([1, 0])
        polymer1.free_energy = -1.0

        polymer2 = Mock(spec=Polymer)
        polymer2.monomer_counts = np.array([0, 1])
        polymer2.free_energy = -2.0

        polymers = [polymer1, polymer2]

        # Mock the subprocess call
        mock_result = Mock()
        mock_result.returncode = 0

        # Mock output file content
        expected_concentrations = np.array([1.5e-7, 2.3e-8])

        with patch.object(runner, "check_coffee_available", return_value=True), patch(
            "tbnexplorer2.coffee.subprocess.run", return_value=mock_result
        ) as mock_run, patch.object(runner, "_parse_coffee_output", return_value=expected_concentrations):
            result = runner.compute_equilibrium_concentrations(polymers, tbn)

            # Verify subprocess was called WITHOUT temperature parameter (default 37°C)
            assert mock_run.called
            call_args = mock_run.call_args[0][0]  # Get the command list
            assert "--temp" not in call_args  # Should not include --temp for default

            # Verify result
            np.testing.assert_array_almost_equal(result, expected_concentrations)

    def test_compute_equilibrium_concentrations_coffee_error(self):
        """Test error handling when COFFEE fails."""
        runner = COFFEERunner("/path/to/coffee")

        tbn = Mock(spec=TBN)
        tbn.concentrations = np.array([1e-7])  # Already in Molar (equivalent to 100 nM)
        tbn.concentration_units = "nM"

        polymer = Mock(spec=Polymer)
        polymer.monomer_counts = np.array([1])
        polymer.free_energy = -1.0
        polymers = [polymer]

        # Mock subprocess to return error
        mock_result = Mock()
        mock_result.returncode = 1
        mock_result.stderr = "COFFEE error message"

        with patch.object(runner, "check_coffee_available", return_value=True), patch(
            "tbnexplorer2.coffee.subprocess.run", return_value=mock_result
        ):
            with pytest.raises(RuntimeError, match="COFFEE failed"):
                runner.compute_equilibrium_concentrations(polymers, tbn)

    def test_compute_equilibrium_with_temperature(self):
        """Test that non-default temperature parameter is passed to coffee-cli."""
        runner = COFFEERunner("/path/to/coffee", temperature=25.0)  # Non-default temperature

        # Mock TBN
        tbn = Mock(spec=TBN)
        tbn.concentrations = np.array([1e-7])
        tbn.concentration_units = "nM"

        # Mock polymer
        polymer = Mock(spec=Polymer)
        polymer.monomer_counts = np.array([1])
        polymer.compute_free_energy = Mock(return_value=-1.0)
        polymers = [polymer]

        # Mock successful subprocess result
        mock_result = Mock()
        mock_result.returncode = 0

        expected_concentrations = np.array([1.5e-7])

        with patch.object(runner, "check_coffee_available", return_value=True), patch(
            "tbnexplorer2.coffee.subprocess.run", return_value=mock_result
        ) as mock_run, patch.object(runner, "_parse_coffee_output", return_value=expected_concentrations):
            result = runner.compute_equilibrium_concentrations(polymers, tbn)

            # Verify subprocess was called with temperature parameter
            assert mock_run.called
            call_args = mock_run.call_args[0][0]  # Get the command list
            assert "--temp" in call_args
            temp_index = call_args.index("--temp")
            assert call_args[temp_index + 1] == "25.0"  # Non-default temperature

            # Verify result
            np.testing.assert_array_almost_equal(result, expected_concentrations)

    def test_compute_equilibrium_default_temp_no_param(self):
        """Test that default temperature (37°C) does NOT pass --temp parameter for backward compatibility."""
        runner = COFFEERunner("/path/to/coffee")  # Uses default 37°C
        assert runner.temperature == 37.0

        # Mock TBN
        tbn = Mock(spec=TBN)
        tbn.concentrations = np.array([1e-7])
        tbn.concentration_units = "nM"

        # Mock polymer
        polymer = Mock(spec=Polymer)
        polymer.monomer_counts = np.array([1])
        polymer.compute_free_energy = Mock(return_value=-1.0)
        polymers = [polymer]

        # Mock successful subprocess result
        mock_result = Mock()
        mock_result.returncode = 0

        expected_concentrations = np.array([1.5e-7])

        with patch.object(runner, "check_coffee_available", return_value=True), patch(
            "tbnexplorer2.coffee.subprocess.run", return_value=mock_result
        ) as mock_run, patch.object(runner, "_parse_coffee_output", return_value=expected_concentrations):
            result = runner.compute_equilibrium_concentrations(polymers, tbn)

            # Verify subprocess was called WITHOUT --temp parameter (backward compatibility)
            assert mock_run.called
            call_args = mock_run.call_args[0][0]  # Get the command list
            assert "--temp" not in call_args  # Should not include --temp for default 37°C

            # Verify result
            np.testing.assert_array_almost_equal(result, expected_concentrations)
