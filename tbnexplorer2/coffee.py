"""
Integration with COFFEE (Computation Of Free-Energy Equilibria) tool.
"""

import os
import subprocess
import tempfile
from typing import TYPE_CHECKING, List, Optional

import numpy as np

from .config import COFFEE_CLI_PATH
from .model import TBN

if TYPE_CHECKING:
    from .polymer_basis import Polymer


class COFFEERunner:
    """Handles integration with the COFFEE equilibrium solver."""

    def __init__(self, coffee_path: str = COFFEE_CLI_PATH):
        """
        Initialize COFFEE runner.

        Args:
            coffee_path: Path to coffee-cli executable
        """
        self.coffee_path = coffee_path

    def check_coffee_available(self) -> bool:
        """Check if COFFEE executable is available."""
        return os.path.isfile(self.coffee_path) and os.access(self.coffee_path, os.X_OK)

    def compute_equilibrium_concentrations(
        self, polymers: List["Polymer"], tbn: TBN, output_dir: Optional[str] = None
    ) -> np.ndarray:
        """
        Compute equilibrium concentrations for polymers using COFFEE.

        Args:
            polymers: List of Polymer objects
            tbn: TBN model with monomer concentrations
            output_dir: Optional directory for temporary files

        Returns:
            Array of polymer concentrations in same order as input

        Raises:
            ValueError: If TBN doesn't have concentrations
            RuntimeError: If COFFEE computation fails
        """
        if tbn.concentrations is None:
            raise ValueError("Cannot compute equilibrium concentrations without monomer concentrations")

        if not self.check_coffee_available():
            raise RuntimeError(f"COFFEE not found at {self.coffee_path}")

        # Create temporary directory if needed
        use_temp_dir = output_dir is None
        if use_temp_dir:
            temp_dir_obj = tempfile.TemporaryDirectory()
            work_dir = temp_dir_obj.name
        else:
            work_dir = output_dir

        try:
            # Prepare CFE file (polymer matrix with free energies)
            cfe_path = os.path.join(work_dir, "polymers.cfe")
            self._write_cfe_file(polymers, cfe_path)

            # Prepare CON file (monomer concentrations)
            con_path = os.path.join(work_dir, "monomers.con")
            self._write_con_file(tbn, con_path)

            # Run COFFEE
            output_path = os.path.join(work_dir, "equilibrium.txt")
            result = subprocess.run(
                [self.coffee_path, cfe_path, con_path, "-o", output_path], capture_output=True, text=True, check=False
            )

            if result.returncode != 0:
                raise RuntimeError(f"COFFEE failed: {result.stderr}")

            # Parse output
            concentrations = self._parse_coffee_output(output_path)

            if len(concentrations) != len(polymers):
                raise RuntimeError(
                    f"COFFEE output has {len(concentrations)} concentrations but expected {len(polymers)}"
                )

            return concentrations

        finally:
            if use_temp_dir:
                temp_dir_obj.cleanup()

    def _write_cfe_file(self, polymers: List["Polymer"], filepath: str):
        """
        Write CFE file for COFFEE.

        Format: Each line has monomer counts followed by free energy.
        """
        with open(filepath, "w") as f:
            for polymer in polymers:
                # Write monomer counts
                counts_str = " ".join(str(int(c)) for c in polymer.monomer_counts)
                # Compute and write free energy
                free_energy = polymer.compute_free_energy()
                f.write(f"{counts_str} {free_energy}\n")

    def _write_con_file(self, tbn: TBN, filepath: str):
        """
        Write CON file for COFFEE.

        Format: One concentration per line, in order of monomers.
        COFFEE expects concentrations in Molar units.
        """
        with open(filepath, "w") as f:
            for conc in tbn.concentrations:
                # tbn.concentrations already returns values in Molar units
                f.write(f"{conc}\n")

    def _parse_coffee_output(self, filepath: str) -> np.ndarray:
        """
        Parse COFFEE output file.

        Returns:
            Array of concentrations
        """
        with open(filepath) as f:
            content = f.read().strip()

        # Parse space-separated values (may be in scientific notation)
        values = content.split()
        concentrations = []

        for val_str in values:
            try:
                # Handle scientific notation like "4.47e-53" or "0.00e0"
                conc = float(val_str)
                concentrations.append(conc)
            except ValueError as e:
                raise RuntimeError(f"Cannot parse concentration value: {val_str}") from e

        return np.array(concentrations)
