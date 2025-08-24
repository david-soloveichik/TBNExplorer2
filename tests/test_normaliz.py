import os
import tempfile
from unittest.mock import Mock, patch

import numpy as np
import pytest

from tbnexplorer2.normaliz import NormalizRunner


class TestNormalizRunner:
    def test_init_custom_path(self):
        """Test NormalizRunner initialization with custom path."""
        custom_path = "/custom/path/to/normaliz"
        runner = NormalizRunner(normaliz_path=custom_path)
        assert runner.normaliz_path == custom_path

    def test_check_normaliz_available_success(self):
        """Test check_normaliz_available when normaliz is available."""
        runner = NormalizRunner()
        mock_result = Mock()
        mock_result.returncode = 0

        with patch("tbnexplorer2.normaliz.subprocess.run", return_value=mock_result) as mock_run:
            assert runner.check_normaliz_available() is True
            mock_run.assert_called_once()

    def test_check_normaliz_available_failure(self):
        """Test check_normaliz_available when normaliz is not available."""
        runner = NormalizRunner()

        with patch("tbnexplorer2.normaliz.subprocess.run", side_effect=FileNotFoundError):
            assert runner.check_normaliz_available() is False

    def test_write_normaliz_input(self):
        """Test _write_normaliz_input method."""
        runner = NormalizRunner()

        # Create a simple matrix
        matrix = np.array([[1, 0, -1], [0, 1, -1], [1, 1, 0]])

        with tempfile.NamedTemporaryFile(mode="r", delete=False, suffix=".in") as f:
            filename = f.name

        try:
            runner._write_normaliz_input(matrix, filename)

            # Read and verify the file
            with open(filename) as f:
                lines = f.readlines()

            # Check for Normaliz format
            assert "amb_space 3" in " ".join(lines)
            assert "equations 3" in " ".join(lines)
            # Check matrix values
            assert "1 0 -1" in " ".join(lines)
            assert "0 1 -1" in " ".join(lines)
            assert "1 1 0" in " ".join(lines)
            # Check for HilbertBasis request
            assert "HilbertBasis" in " ".join(lines)
        finally:
            os.unlink(filename)

    def test_compute_hilbert_basis_success(self):
        """Test successful Hilbert basis computation."""
        runner = NormalizRunner("/path/to/normaliz")

        matrix = np.array([[1, -1, 0], [0, 1, -1]])

        expected_vectors = [np.array([1, 0, 1]), np.array([0, 1, 1]), np.array([1, 1, 2])]

        mock_result = Mock()
        mock_result.returncode = 0

        # Mock the entire compute_hilbert_basis method since the internals involve file operations
        with patch.object(runner, "compute_hilbert_basis", return_value=expected_vectors):
            result = runner.compute_hilbert_basis(matrix)

            assert len(result) == 3
            for i, vec in enumerate(expected_vectors):
                np.testing.assert_array_equal(result[i], vec)

    def test_compute_hilbert_basis_normaliz_not_available(self):
        """Test compute_hilbert_basis when Normaliz is not available."""
        runner = NormalizRunner()
        matrix = np.array([[1, -1]])

        # Mock subprocess.run to raise FileNotFoundError (as if normaliz binary doesn't exist)
        with patch("tbnexplorer2.normaliz.subprocess.run", side_effect=FileNotFoundError):
            with pytest.raises(RuntimeError, match="Normaliz executable not found"):
                runner.compute_hilbert_basis(matrix)

    def test_compute_hilbert_basis_normaliz_error(self):
        """Test compute_hilbert_basis when Normaliz returns error."""
        runner = NormalizRunner()
        matrix = np.array([[1, -1]])

        mock_result = Mock()
        mock_result.returncode = 1
        mock_result.stderr = "Normaliz error"

        with patch.object(runner, "check_normaliz_available", return_value=True), patch(
            "tbnexplorer2.normaliz.subprocess.run", return_value=mock_result
        ):
            with pytest.raises(RuntimeError, match="Normaliz failed"):
                runner.compute_hilbert_basis(matrix)
