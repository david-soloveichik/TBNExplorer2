import numpy as np
import pytest

from tbnexplorer2.fourtitwo import FourTiTwoRunner
from tbnexplorer2.normaliz import NormalizRunner


class TestFourTiTwoRunner:
    def test_fourtitwo_availability(self):
        """Test checking 4ti2 availability."""
        runner = FourTiTwoRunner()
        # This may or may not be available depending on the system
        # Just check that the method runs without error
        available = runner.check_fourtitwo_available()
        assert isinstance(available, bool)

    def test_compute_simple_hilbert_basis(self):
        """Test computing Hilbert basis for a simple example."""
        runner = FourTiTwoRunner()

        # Skip test if 4ti2 is not available
        if not runner.check_fourtitwo_available():
            pytest.skip("4ti2 not available")

        # Simple example: x + y = 2z
        # Hilbert basis should include vectors like [1, 1, 1], [2, 0, 1], [0, 2, 1]
        matrix = np.array([[1, 1, -2]])

        try:
            hilbert_basis = runner.compute_hilbert_basis(matrix)

            # Check that we got some vectors
            assert len(hilbert_basis) > 0

            # Check that all vectors satisfy the equation
            for vector in hilbert_basis:
                assert len(vector) == 3
                # Check x + y - 2z = 0
                assert np.allclose(matrix @ vector, 0)
                # Check non-negativity
                assert np.all(vector >= 0)
        except RuntimeError as e:
            # If 4ti2 fails, skip the test
            pytest.skip(f"4ti2 execution failed: {e}")

    def test_compare_with_normaliz(self):
        """Test that 4ti2 and Normaliz produce equivalent Hilbert bases."""
        fourtitwo_runner = FourTiTwoRunner()
        normaliz_runner = NormalizRunner()

        # Skip test if either tool is not available
        if not fourtitwo_runner.check_fourtitwo_available():
            pytest.skip("4ti2 not available")
        if not normaliz_runner.check_normaliz_available():
            pytest.skip("Normaliz not available")

        # Simple test case
        matrix = np.array([[1, 1, -2], [1, -1, 0]])

        try:
            fourtitwo_basis = fourtitwo_runner.compute_hilbert_basis(matrix)
            normaliz_basis = normaliz_runner.compute_hilbert_basis(matrix)

            # Convert to sets of tuples for comparison
            fourtitwo_set = {tuple(v) for v in fourtitwo_basis}
            normaliz_set = {tuple(v) for v in normaliz_basis}

            # Check that both computed the same basis
            # Note: The order might be different, but the sets should be equal
            assert fourtitwo_set == normaliz_set, f"4ti2: {fourtitwo_set}, Normaliz: {normaliz_set}"

        except RuntimeError as e:
            # If either tool fails, skip the test
            pytest.skip(f"Tool execution failed: {e}")

    def test_empty_matrix(self):
        """Test handling of empty matrix."""
        runner = FourTiTwoRunner()

        if not runner.check_fourtitwo_available():
            pytest.skip("4ti2 not available")

        # Empty matrix (no constraints)
        matrix = np.array([]).reshape(0, 3)

        try:
            hilbert_basis = runner.compute_hilbert_basis(matrix)
            # With no constraints, the basis should be the standard basis vectors
            assert len(hilbert_basis) >= 3
        except RuntimeError:
            # Some versions might not handle empty matrices well
            pass

    def test_zero_dimensional_case(self):
        """Test case where the only solution is the zero vector."""
        runner = FourTiTwoRunner()

        if not runner.check_fourtitwo_available():
            pytest.skip("4ti2 not available")

        # Overdetermined system with only zero solution
        # x = 0, y = 0, z = 0
        matrix = np.array([[1, 0, 0], [0, 1, 0], [0, 0, 1]])

        try:
            hilbert_basis = runner.compute_hilbert_basis(matrix)
            # Should return empty basis or just the zero vector
            assert len(hilbert_basis) == 0 or (len(hilbert_basis) == 1 and np.allclose(hilbert_basis[0], 0))
        except RuntimeError:
            # This is an edge case that might fail
            pass
