import os
import subprocess
import tempfile
from typing import List

import numpy as np

from .config import FOURTI2_PATH


class FourTiTwoRunner:
    """Wrapper for running 4ti2 and parsing its output."""

    def __init__(self, fourtitwo_path: str = FOURTI2_PATH):
        """
        Initialize 4ti2 runner.

        Args:
            fourtitwo_path: Path to 4ti2 installation directory
        """
        self.fourtitwo_path = fourtitwo_path
        self.hilbert_executable = os.path.join(fourtitwo_path, "bin", "hilbert")
        self.zsolve_executable = os.path.join(fourtitwo_path, "bin", "zsolve")

    def compute_hilbert_basis(self, matrix: np.ndarray) -> List[np.ndarray]:
        """
        Compute Hilbert basis of the cone {x >= 0 : matrix * x = 0}.

        Args:
            matrix: Matrix defining the linear equations

        Returns:
            List of Hilbert basis vectors

        Raises:
            RuntimeError: If 4ti2 execution fails
        """
        # Create temporary directory for 4ti2 files
        with tempfile.TemporaryDirectory() as tmpdir:
            base_name = os.path.join(tmpdir, "problem")

            # Write 4ti2 input files
            self._write_fourtitwo_input(matrix, base_name)

            # Try hilbert first, fall back to zsolve if needed
            try:
                # Run hilbert
                output_file = self._run_hilbert(base_name)
                # Parse Hilbert basis from output
                hilbert_basis = self._parse_hilbert_output(output_file)
            except (RuntimeError, FileNotFoundError):
                # If hilbert fails, try zsolve
                output_file = self._run_zsolve(base_name)
                # Parse Hilbert basis from zsolve output
                hilbert_basis = self._parse_zsolve_output(output_file)

        return hilbert_basis

    def _write_fourtitwo_input(self, matrix: np.ndarray, base_name: str):
        """
        Write 4ti2 input files for computing Hilbert basis.

        Args:
            matrix: Matrix defining equations (each row is an equation)
            base_name: Base path for input files (without extension)
        """
        n_equations, n_variables = matrix.shape

        # Write the .mat file (equation matrix)
        mat_file = base_name + ".mat"
        with open(mat_file, "w") as f:
            # First line: dimensions (rows columns)
            f.write(f"{n_equations} {n_variables}\n")
            # Write matrix rows
            for row in matrix:
                f.write(" ".join(str(int(val)) for val in row) + "\n")

        # Write the .sign file (all variables non-negative)
        sign_file = base_name + ".sign"
        with open(sign_file, "w") as f:
            f.write(f"1 {n_variables}\n")
            f.write(" ".join("+" for _ in range(n_variables)) + "\n")

        # Write the .rel file (all equations)
        rel_file = base_name + ".rel"
        with open(rel_file, "w") as f:
            f.write(f"1 {n_equations}\n")
            f.write(" ".join("=" for _ in range(n_equations)) + "\n")

    def _run_hilbert(self, base_name: str) -> str:
        """
        Run 4ti2 hilbert on the input files.

        Args:
            base_name: Base path for input/output files

        Returns:
            Path to output file

        Raises:
            RuntimeError: If hilbert execution fails
        """
        try:
            # Run hilbert
            result = subprocess.run([self.hilbert_executable, base_name], capture_output=True, text=True, check=False)

            if result.returncode != 0:
                raise RuntimeError(f"4ti2 hilbert failed: {result.stderr}")

            # hilbert creates output file with .hil extension
            output_file = base_name + ".hil"

            if not os.path.exists(output_file):
                raise RuntimeError(f"4ti2 hilbert output file not found: {output_file}")

            return output_file

        except FileNotFoundError as e:
            raise RuntimeError(
                f"4ti2 hilbert executable not found at '{self.hilbert_executable}'. "
                f"Please install 4ti2 or update FOURTI2_PATH"
            ) from e

    def _run_zsolve(self, base_name: str) -> str:
        """
        Run 4ti2 zsolve on the input files.

        Args:
            base_name: Base path for input/output files

        Returns:
            Path to output file

        Raises:
            RuntimeError: If zsolve execution fails
        """
        try:
            # Run zsolve
            result = subprocess.run([self.zsolve_executable, base_name], capture_output=True, text=True, check=False)

            if result.returncode != 0:
                raise RuntimeError(f"4ti2 zsolve failed: {result.stderr}")

            # zsolve creates output file with .zhom extension for homogeneous solutions
            output_file = base_name + ".zhom"

            if not os.path.exists(output_file):
                # Try .zinhom for inhomogeneous solutions as fallback
                output_file = base_name + ".zinhom"
                if not os.path.exists(output_file):
                    raise RuntimeError("4ti2 zsolve output file not found")

            return output_file

        except FileNotFoundError as e:
            raise RuntimeError(
                f"4ti2 zsolve executable not found at '{self.zsolve_executable}'. "
                f"Please install 4ti2 or update FOURTI2_PATH"
            ) from e

    def _parse_hilbert_output(self, output_file: str) -> List[np.ndarray]:
        """
        Parse Hilbert basis from 4ti2 hilbert output file.

        Args:
            output_file: Path to hilbert output file (.hil)

        Returns:
            List of Hilbert basis vectors as numpy arrays
        """
        hilbert_basis = []

        with open(output_file) as f:
            lines = f.readlines()

            # First line contains dimensions: n_vectors n_variables
            if not lines:
                return hilbert_basis

            first_line = lines[0].strip().split()
            if len(first_line) != 2:
                raise RuntimeError(f"Invalid 4ti2 hilbert output format: {lines[0]}")

            n_vectors = int(first_line[0])
            n_variables = int(first_line[1])

            # Read the vectors
            for i in range(1, min(n_vectors + 1, len(lines))):
                line = lines[i].strip()
                if line:
                    try:
                        vector = [int(x) for x in line.split()]
                        if len(vector) == n_variables:
                            hilbert_basis.append(np.array(vector))
                    except ValueError:
                        pass  # Skip lines that can't be parsed

        return hilbert_basis

    def _parse_zsolve_output(self, output_file: str) -> List[np.ndarray]:
        """
        Parse Hilbert basis from 4ti2 zsolve output file.

        Args:
            output_file: Path to zsolve output file (.zhom or .zinhom)

        Returns:
            List of Hilbert basis vectors as numpy arrays
        """
        # zsolve output format is similar to hilbert
        return self._parse_hilbert_output(output_file)

    def compute_module_generators_for_slice(self, equations: np.ndarray, slice_vector: np.ndarray) -> List[np.ndarray]:
        """
        Compute module generators over original monoid for a slice using 4ti2's zsolve.

        Computes the minimal inhomogeneous solutions (module generators) of:
        { x >= 0 : equations * x = 0, slice_vector * x >= 1 }

        This is used for finding reactions that produce a target polymer,
        where slice_vector has 1 at the target polymer position and 0 elsewhere.

        Args:
            equations: Matrix defining linear equations (B matrix for mass conservation)
            slice_vector: Row vector defining the slice (e.g., selecting a target polymer)

        Returns:
            List of module generator vectors

        Raises:
            RuntimeError: If 4ti2 execution fails
        """
        # Create temporary directory for 4ti2 files
        with tempfile.TemporaryDirectory() as tmpdir:
            base_name = os.path.join(tmpdir, "slice")

            # Write 4ti2 input files for the slice problem
            self._write_zsolve_slice_input(equations, slice_vector, base_name)

            # Run zsolve
            try:
                result = subprocess.run(
                    [self.zsolve_executable, base_name], capture_output=True, text=True, check=False
                )

                if result.returncode != 0:
                    raise RuntimeError(f"4ti2 zsolve failed: {result.stderr}")

                # Parse the inhomogeneous solutions (.zinhom file)
                output_file = base_name + ".zinhom"

                if not os.path.exists(output_file):
                    # No inhomogeneous solutions found (empty result)
                    return []

                # Parse module generators from output
                module_generators = self._parse_zsolve_output(output_file)

            except FileNotFoundError as e:
                raise RuntimeError(
                    f"4ti2 zsolve executable not found at '{self.zsolve_executable}'. "
                    f"Please install 4ti2 or update FOURTI2_PATH"
                ) from e

        return module_generators

    def _write_zsolve_slice_input(self, equations: np.ndarray, slice_vector: np.ndarray, base_name: str):
        """
        Write 4ti2 zsolve input files for computing module generators for a slice.

        We solve: { x >= 0 : equations * x = 0, slice * x >= 1 }

        According to the 4ti2 documentation format:
        - .mat file: coefficient matrix (equations followed by slice row)
        - .rel file: relation types ('=' for equations, '>' for slice inequality)
        - .rhs file: right-hand sides (0s for equations, 1 for slice)
        - .sign file: variable sign restrictions (0 for free variables)

        Args:
            equations: Matrix of equations (each row is an equation)
            slice_vector: Row vector defining the slice
            base_name: Base path for input files (without extension)
        """
        n_equations, n_variables = equations.shape

        # Total number of rows: equations + slice constraint
        n_rows = n_equations + 1

        # Write the .mat file (coefficient matrix)
        mat_file = base_name + ".mat"
        with open(mat_file, "w") as f:
            f.write(f"{n_rows} {n_variables}\n")
            # Write equation rows
            for row in equations:
                f.write(" ".join(str(int(val)) for val in row) + "\n")
            # Write slice row
            f.write(" ".join(str(int(val)) for val in slice_vector) + "\n")

        # Write the .rel file (relation types)
        rel_file = base_name + ".rel"
        with open(rel_file, "w") as f:
            f.write(f"1 {n_rows}\n")
            # '=' for equations, '>' for slice inequality
            relations = ["="] * n_equations + [">"]
            f.write(" ".join(relations) + "\n")

        # Write the .rhs file (right-hand sides)
        rhs_file = base_name + ".rhs"
        with open(rhs_file, "w") as f:
            f.write(f"1 {n_rows}\n")
            # 0 for equations, 1 for slice constraint
            rhs_values = ["0"] * n_equations + ["1"]
            f.write(" ".join(rhs_values) + "\n")

        # Write the .sign file (variable sign restrictions)
        # Use 1 for non-negative variables (since we're in lifted space)
        # This ensures all variables are >= 0
        sign_file = base_name + ".sign"
        with open(sign_file, "w") as f:
            f.write(f"1 {n_variables}\n")
            f.write(" ".join("1" for _ in range(n_variables)) + "\n")

    def check_fourtitwo_available(self) -> bool:
        """
        Check if 4ti2 is available at the configured path.

        Returns:
            True if 4ti2 is available, False otherwise
        """
        # Check if either hilbert or zsolve is available
        hilbert_available = os.path.exists(self.hilbert_executable) and os.access(self.hilbert_executable, os.X_OK)
        zsolve_available = os.path.exists(self.zsolve_executable) and os.access(self.zsolve_executable, os.X_OK)
        return hilbert_available or zsolve_available
