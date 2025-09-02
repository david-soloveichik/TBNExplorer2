"""
Integration with NUPACK concentrations tool.
"""

import os
import subprocess
import tempfile
from typing import TYPE_CHECKING, List, Optional

import numpy as np

from .config import NUPACK_CONCENTRATIONS_PATH
from .model import TBN

if TYPE_CHECKING:
    from .polymer_basis import Polymer


class NupackRunner:
    """Handles integration with the NUPACK concentrations equilibrium solver."""

    def __init__(self, nupack_path: str = NUPACK_CONCENTRATIONS_PATH, temperature: float = 37.0):
        """
        Initialize NUPACK runner.

        Args:
            nupack_path: Path to NUPACK concentrations executable
            temperature: Temperature in Celsius (defaults to 37)
        """
        self.nupack_path = nupack_path
        self.temperature = temperature

    def check_nupack_available(self) -> bool:
        """Check if NUPACK executable is available."""
        return os.path.isfile(self.nupack_path) and os.access(self.nupack_path, os.X_OK)

    def compute_equilibrium_concentrations(
        self,
        polymers: List["Polymer"],
        tbn: TBN,
        output_dir: Optional[str] = None,
        deltaG: Optional[List[float]] = None,
        temperature: float = 37.0,
    ) -> np.ndarray:
        """
        Compute equilibrium concentrations for polymers using NUPACK.

        Args:
            polymers: List of Polymer objects
            tbn: TBN model with monomer concentrations
            output_dir: Optional directory for temporary files
            deltaG: List of [dG_bond, dG_assoc, dH_assoc] (default: [-1.0, 0.0, 0.0])
            temperature: Temperature in Celsius (default: 37.0)

        Returns:
            Array of polymer concentrations in same order as input

        Raises:
            ValueError: If TBN doesn't have concentrations
            RuntimeError: If NUPACK computation fails
        """
        if tbn.concentrations is None:
            raise ValueError("Cannot compute equilibrium concentrations without monomer concentrations")

        if not self.check_nupack_available():
            raise RuntimeError(f"NUPACK not found at {self.nupack_path}")

        # Create temporary directory if needed
        use_temp_dir = output_dir is None
        if use_temp_dir:
            temp_dir_obj = tempfile.TemporaryDirectory()
            work_dir = temp_dir_obj.name
        else:
            work_dir = output_dir

        try:
            base_name = "nupack_input"
            base_path = os.path.join(work_dir, base_name)

            # Prepare OCX file (polymer matrix with free energies)
            ocx_path = f"{base_path}.ocx"
            self._write_ocx_file(polymers, ocx_path, deltaG, temperature)

            # Prepare CON file (monomer concentrations)
            con_path = f"{base_path}.con"
            self._write_con_file(tbn, con_path)

            # Run NUPACK concentrations
            # Command: concentrations -sort 0 -T {temperature_in_C} {base}
            # -sort 0 preserves input order (no sorting)
            cmd = [self.nupack_path, "-sort", "0", "-T", str(temperature), base_path]

            result = subprocess.run(cmd, capture_output=True, text=True, check=False, cwd=work_dir)

            if result.returncode != 0:
                raise RuntimeError(f"NUPACK failed: {result.stderr}")

            # Parse output from .eq file
            eq_path = f"{base_path}.eq"
            concentrations = self._parse_nupack_output(eq_path)

            if len(concentrations) != len(polymers):
                raise RuntimeError(
                    f"NUPACK output has {len(concentrations)} concentrations but expected {len(polymers)}"
                )

            return concentrations

        finally:
            if use_temp_dir:
                temp_dir_obj.cleanup()

    def _write_ocx_file(
        self, polymers: List["Polymer"], filepath: str, deltaG: Optional[List[float]] = None, temperature: float = 37.0
    ):
        """
        Write OCX file for NUPACK.

        Format: Tab-delimited file with columns:
        - polymer id (line number starting from 1)
        - always 1
        - monomer counts
        - free energy

        Args:
            polymers: List of Polymer objects
            filepath: Path to write OCX file
            deltaG: List of [dG_bond, dG_assoc, dH_assoc] (default: [-1.0, 0.0, 0.0])
            temperature: Temperature in Celsius (default: 37.0)
        """
        with open(filepath, "w") as f:
            for idx, polymer in enumerate(polymers, 1):
                # Polymer ID (line number)
                row = [str(idx)]
                # Always 1
                row.append("1")
                # Monomer counts
                row.extend(str(int(c)) for c in polymer.monomer_counts)
                # Free energy
                free_energy = polymer.compute_free_energy(deltaG, temperature)
                row.append(str(free_energy))
                # Write tab-delimited row
                f.write("\t".join(row) + "\n")

    def _write_con_file(self, tbn: TBN, filepath: str):
        """
        Write CON file for NUPACK.

        Format: One concentration per line in Molar units, in order of monomers.
        """
        with open(filepath, "w") as f:
            for conc in tbn.concentrations:
                # tbn.concentrations already returns values in Molar units
                f.write(f"{conc}\n")

    def _parse_nupack_output(self, filepath: str) -> np.ndarray:
        """
        Parse NUPACK .eq output file.

        The .eq file adds a concentration column (last column) to the .ocx format.

        Returns:
            Array of concentrations in Molar units
        """
        concentrations = []

        with open(filepath) as f:
            for line in f:
                line = line.strip()
                if not line:
                    continue
                # Skip header lines (e.g., "% NUPACK 3.2.2")
                if line.startswith("%") or line.startswith("#"):
                    continue
                # Split by tabs and get the last column (concentration)
                parts = line.split("\t")
                if len(parts) < 2:
                    # Skip lines that don't have enough columns
                    continue
                # Last column is the concentration
                conc_str = parts[-1]
                try:
                    conc = float(conc_str)
                    concentrations.append(conc)
                except ValueError:
                    # Skip lines where last column is not a number
                    continue

        arr = np.array(concentrations)
        if arr.size == 0:
            raise RuntimeError("Invalid .eq file format: no concentration data parsed")
        return arr
