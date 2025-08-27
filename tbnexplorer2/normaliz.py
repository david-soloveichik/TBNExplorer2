import os
import subprocess
import tempfile
from typing import List

import numpy as np

from .config import NORMALIZ_PATH


class NormalizRunner:
    """Wrapper for running Normaliz and parsing its output."""

    def __init__(self, normaliz_path: str = NORMALIZ_PATH):
        """
        Initialize Normaliz runner.

        Args:
            normaliz_path: Path to Normaliz executable
        """
        self.normaliz_path = normaliz_path

    def compute_hilbert_basis(self, matrix: np.ndarray) -> List[np.ndarray]:
        """
        Compute Hilbert basis of the cone {x >= 0 : matrix * x = 0}.

        Args:
            matrix: Matrix defining the linear equations

        Returns:
            List of Hilbert basis vectors

        Raises:
            RuntimeError: If Normaliz execution fails
        """
        # Create temporary directory for Normaliz files
        with tempfile.TemporaryDirectory() as tmpdir:
            input_file = os.path.join(tmpdir, "input.in")

            # Write Normaliz input file
            self._write_normaliz_input(matrix, input_file)

            # Run Normaliz
            output_file = self._run_normaliz(input_file)

            # Parse Hilbert basis from output
            hilbert_basis = self._parse_hilbert_basis(output_file)

        return hilbert_basis

    def _write_normaliz_input(self, matrix: np.ndarray, filepath: str):
        """
        Write Normaliz input file for computing Hilbert basis.

        Args:
            matrix: Matrix defining equations (each row is an equation)
            filepath: Path to write input file
        """
        n_equations, n_variables = matrix.shape

        with open(filepath, "w") as f:
            f.write("/* Normaliz input for Hilbert basis computation */\n\n")

            # Specify ambient space dimension
            f.write(f"amb_space {n_variables}\n\n")

            # Write equations
            f.write(f"equations {n_equations}\n")
            for row in matrix:
                f.write(" ".join(str(int(val)) for val in row) + "\n")
            f.write("\n")

            # Request Hilbert basis computation
            f.write("HilbertBasis\n")

    def _run_normaliz(self, input_file: str) -> str:
        """
        Run Normaliz on the input file.

        Args:
            input_file: Path to Normaliz input file

        Returns:
            Path to output file

        Raises:
            RuntimeError: If Normaliz execution fails
        """
        try:
            # Run Normaliz
            result = subprocess.run([self.normaliz_path, input_file], capture_output=True, text=True, check=False)

            if result.returncode != 0:
                error_msg = result.stderr if result.stderr else result.stdout
                raise RuntimeError(f"Normaliz failed with return code {result.returncode}: {error_msg}")

            # Normaliz creates output file with .out extension
            output_file = input_file.rsplit(".", 1)[0] + ".out"

            if not os.path.exists(output_file):
                raise RuntimeError(f"Normaliz output file not found: {output_file}")

            return output_file

        except FileNotFoundError as e:
            raise RuntimeError(
                f"Normaliz executable not found at '{self.normaliz_path}'. "
                f"Please install Normaliz or update NORMALIZ_PATH in normaliz.py"
            ) from e

    def compute_module_generators_for_slice(self, equations: np.ndarray, slice_vector: np.ndarray) -> List[np.ndarray]:
        """
        Compute module generators over original monoid for a slice.

        Computes the module generators of the slice:
        { x >= 0 : equations * x = 0, slice_vector * x >= 1 }

        This is specifically for finding reactions that produce a target polymer,
        where slice_vector has 1 at the target polymer position and 0 elsewhere.

        Args:
            equations: Matrix defining linear equations (B matrix for mass conservation)
            slice_vector: Row vector defining the slice (e.g., selecting a target polymer)

        Returns:
            List of module generator vectors

        Raises:
            RuntimeError: If Normaliz execution fails
        """
        # Create temporary directory for Normaliz files
        with tempfile.TemporaryDirectory() as tmpdir:
            input_file = os.path.join(tmpdir, "input.in")

            # Write Normaliz input file for module generators
            self._write_normaliz_input_for_module_generators(equations, slice_vector, input_file)

            # Run Normaliz
            output_file = self._run_normaliz(input_file)

            # Parse module generators from output
            module_generators = self._parse_module_generators(output_file)

        return module_generators

    def _write_normaliz_input_for_module_generators(
        self, equations: np.ndarray, slice_vector: np.ndarray, filepath: str
    ):
        """
        Write Normaliz input file for computing module generators for a slice.

        We solve: { x >= 0 : equations * x = 0, slice * x >= 1 }
        This represents reactions that conserve mass and produce at least one target polymer.

        Args:
            equations: Matrix of equations (each row is an equation)
            slice_vector: Row vector defining the slice (positivity constraint)
            filepath: Path to write input file
        """
        n_variables = equations.shape[1]

        with open(filepath, "w") as f:
            f.write("/* Normaliz input for Hilbert basis with strict inequality */\n\n")

            # Specify ambient space dimension
            f.write(f"amb_space {n_variables}\n\n")

            # Write equations if any
            if equations.shape[0] > 0:
                f.write(f"equations {equations.shape[0]}\n")
                for row in equations:
                    f.write(" ".join(str(int(val)) for val in row) + "\n")
                f.write("\n")

            # Use strict_inequalities to enforce slice * x >= 1
            # This ensures the target polymer must be produced
            f.write("strict_inequalities 1\n")
            f.write(" ".join(str(int(val)) for val in slice_vector) + "\n")
            f.write("\n")

            # Request Hilbert basis computation for the slice
            # This gives us all solutions, not just module generators
            f.write("HilbertBasis\n")

    def _parse_module_generators(self, output_file: str) -> List[np.ndarray]:
        """
        Parse module generators (Hilbert basis) from Normaliz output file.

        Args:
            output_file: Path to Normaliz output file

        Returns:
            List of module generator vectors as numpy arrays
        """
        # Use the existing Hilbert basis parser since we're now computing Hilbert basis
        return self._parse_hilbert_basis(output_file)

    def _write_normaliz_input_with_strict(
        self, equations: np.ndarray, inequalities: np.ndarray, strict_inequality: np.ndarray, filepath: str
    ):
        """
        Write Normaliz input file with equations, inequalities, and strict inequality.

        Args:
            equations: Matrix of equations (each row is an equation)
            inequalities: Matrix of inequalities (each row is an inequality)
            strict_inequality: Row vector for strict inequality
            filepath: Path to write input file
        """
        n_variables = equations.shape[1]

        with open(filepath, "w") as f:
            f.write("/* Normaliz input with strict inequality */\n\n")

            # Specify ambient space dimension
            f.write(f"amb_space {n_variables}\n\n")

            # Write equations if any
            if equations.shape[0] > 0:
                f.write(f"equations {equations.shape[0]}\n")
                for row in equations:
                    f.write(" ".join(str(int(val)) for val in row) + "\n")
                f.write("\n")

            # Write non-strict inequalities if any
            if inequalities.shape[0] > 0:
                f.write(f"inequalities {inequalities.shape[0]}\n")
                for row in inequalities:
                    f.write(" ".join(str(int(val)) for val in row) + "\n")
                f.write("\n")

            # Write strict inequality as P*x >= 1
            # Normaliz's strict_inequalities type encodes ξ·x >= 1
            f.write("strict_inequalities 1\n")
            f.write(" ".join(str(int(val)) for val in strict_inequality) + "\n")
            f.write("\n")

            # Request Hilbert basis computation
            f.write("HilbertBasis\n")

    def _parse_hilbert_basis(self, output_file: str) -> List[np.ndarray]:
        """
        Parse Hilbert basis from Normaliz output file.

        This method handles various Normaliz output formats based on the
        implementation in TBNCanonicalReactionsEnumerator.

        Args:
            output_file: Path to Normaliz output file

        Returns:
            List of Hilbert basis vectors as numpy arrays
        """
        hilbert_basis = []
        in_hilbert_section = False
        found_header = False

        with open(output_file) as f:
            for line in f:
                line = line.strip()

                # Look for Hilbert basis section - multiple possible formats
                if (
                    "lattice points in polytope (Hilbert basis elements of degree 1):" in line
                    or "Hilbert basis elements:" in line
                    or "module generators:" in line
                ):
                    in_hilbert_section = True
                    found_header = True
                    continue
                elif "Hilbert basis elements of higher degree:" in line and found_header:
                    # Additional basis elements section
                    continue
                elif in_hilbert_section and any(
                    marker in line
                    for marker in [
                        "extreme rays:",
                        "support hyperplanes:",
                        "equations:",
                        "basis elements of generated",
                        "***",
                    ]
                ):
                    # End of Hilbert basis section
                    break

                # Parse Hilbert basis vectors
                if in_hilbert_section and line and not line.startswith("*") and all(c in "0123456789 -" for c in line):
                    try:
                        vector = [int(x) for x in line.split()]
                        if vector:  # Non-empty vector
                            hilbert_basis.append(np.array(vector))
                    except ValueError:
                        pass  # Skip lines that can't be parsed

        return hilbert_basis

    def check_normaliz_available(self) -> bool:
        """
        Check if Normaliz is available at the configured path.

        Returns:
            True if Normaliz is available, False otherwise
        """
        try:
            result = subprocess.run([self.normaliz_path, "--version"], capture_output=True, text=True, timeout=5)
            return result.returncode == 0
        except (FileNotFoundError, subprocess.TimeoutExpired):
            return False
