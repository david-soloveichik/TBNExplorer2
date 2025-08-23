import os
import tempfile
from unittest.mock import Mock, patch

import numpy as np
import pytest

from tbnexplorer2.normaliz import NormalizRunner


class TestNormalizRunner:
    def test_init_default_path(self):
        """Test NormalizRunner initialization with default path."""
        runner = NormalizRunner()
        assert runner.normaliz_path == "normaliz"  # Default from config

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

        with patch("subprocess.run", return_value=mock_result) as mock_run:
            assert runner.check_normaliz_available() is True
            mock_run.assert_called_once()

    def test_check_normaliz_available_failure(self):
        """Test check_normaliz_available when normaliz is not available."""
        runner = NormalizRunner()

        with patch("subprocess.run", side_effect=FileNotFoundError):
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

            # Check matrix dimensions
            assert "3 3" in lines[0]
            # Check matrix values
            assert "1 0 -1" in lines[1]
            assert "0 1 -1" in lines[2]
            assert "1 1 0" in lines[3]
            # Check equations specification
            assert "equations" in lines[4].lower()
        finally:
            os.unlink(filename)

    def test_parse_normaliz_output_with_vectors(self):
        """Test _parse_normaliz_output with valid vectors."""
        runner = NormalizRunner()

        # Create sample Normaliz output
        output_content = """3 Hilbert basis elements:
 1 0 1
 0 1 1
 1 1 0
        
3 Hilbert basis elements of degree 1:
 1 0 1
 0 1 1
 1 1 0"""

        with tempfile.NamedTemporaryFile(mode="w", delete=False, suffix=".out") as f:
            f.write(output_content)
            filename = f.name

        try:
            vectors = runner._parse_normaliz_output(filename, 3)

            assert len(vectors) == 3
            np.testing.assert_array_equal(vectors[0], [1, 0, 1])
            np.testing.assert_array_equal(vectors[1], [0, 1, 1])
            np.testing.assert_array_equal(vectors[2], [1, 1, 0])
        finally:
            os.unlink(filename)

    def test_parse_normaliz_output_no_hilbert_basis(self):
        """Test _parse_normaliz_output when no Hilbert basis found."""
        runner = NormalizRunner()

        output_content = """Some other output
No Hilbert basis here"""

        with tempfile.NamedTemporaryFile(mode="w", delete=False, suffix=".out") as f:
            f.write(output_content)
            filename = f.name

        try:
            vectors = runner._parse_normaliz_output(filename, 2)
            assert len(vectors) == 0
        finally:
            os.unlink(filename)

    def test_compute_hilbert_basis_success(self):
        """Test successful Hilbert basis computation."""
        runner = NormalizRunner("/path/to/normaliz")

        matrix = np.array([[1, -1, 0], [0, 1, -1]])

        expected_vectors = [np.array([1, 0, 1]), np.array([0, 1, 1]), np.array([1, 1, 2])]

        mock_result = Mock()
        mock_result.returncode = 0

        with patch.object(runner, "check_normaliz_available", return_value=True), patch(
            "subprocess.run", return_value=mock_result
        ), patch.object(runner, "_parse_normaliz_output", return_value=expected_vectors):
            result = runner.compute_hilbert_basis(matrix)

            assert len(result) == 3
            for i, vec in enumerate(expected_vectors):
                np.testing.assert_array_equal(result[i], vec)

    def test_compute_hilbert_basis_normaliz_not_available(self):
        """Test compute_hilbert_basis when Normaliz is not available."""
        runner = NormalizRunner()
        matrix = np.array([[1, -1]])

        with patch.object(runner, "check_normaliz_available", return_value=False):
            with pytest.raises(RuntimeError, match="Normaliz not found"):
                runner.compute_hilbert_basis(matrix)

    def test_compute_hilbert_basis_normaliz_error(self):
        """Test compute_hilbert_basis when Normaliz returns error."""
        runner = NormalizRunner()
        matrix = np.array([[1, -1]])

        mock_result = Mock()
        mock_result.returncode = 1
        mock_result.stderr = "Normaliz error"

        with patch.object(runner, "check_normaliz_available", return_value=True), patch(
            "subprocess.run", return_value=mock_result
        ):
            with pytest.raises(RuntimeError, match="Normaliz computation failed"):
                runner.compute_hilbert_basis(matrix)

    def test_compute_hilbert_basis_with_output_dir(self):
        """Test compute_hilbert_basis with specified output directory."""
        runner = NormalizRunner()
        matrix = np.array([[1, -1, 0]])

        with tempfile.TemporaryDirectory() as temp_dir:
            mock_result = Mock()
            mock_result.returncode = 0

            expected_vectors = [np.array([1, 0, 1])]

            with patch.object(runner, "check_normaliz_available", return_value=True), patch(
                "subprocess.run", return_value=mock_result
            ) as mock_run, patch.object(runner, "_parse_normaliz_output", return_value=expected_vectors):
                result = runner.compute_hilbert_basis(matrix, output_dir=temp_dir)

                # Verify subprocess was called with correct working directory
                call_kwargs = mock_run.call_args[1]
                assert call_kwargs.get("cwd") == temp_dir

                assert len(result) == 1
                np.testing.assert_array_equal(result[0], expected_vectors[0])

    def test_parse_normaliz_output_with_empty_lines(self):
        """Test _parse_normaliz_output handles empty lines correctly."""
        runner = NormalizRunner()

        output_content = """2 Hilbert basis elements:

 1 0
 
 0 1

"""

        with tempfile.NamedTemporaryFile(mode="w", delete=False, suffix=".out") as f:
            f.write(output_content)
            filename = f.name

        try:
            vectors = runner._parse_normaliz_output(filename, 2)

            assert len(vectors) == 2
            np.testing.assert_array_equal(vectors[0], [1, 0])
            np.testing.assert_array_equal(vectors[1], [0, 1])
        finally:
            os.unlink(filename)
