import math
import os
from typing import List, Optional, Tuple

import numpy as np

from .coffee import COFFEERunner
from .model import TBN, Monomer
from .normaliz import NormalizRunner
from .polymat_io import PolymatData, PolymatWriter, check_matrix_hash
from .tbnpolys_io import TbnpolysWriter
from .units import from_molar, get_unit_display_name

# Boltzmann constant in kcal/mol/K
KB = 0.001987204259


def _celcius_to_kelvin(temp_c: float) -> float:
    """Convert temperature from Celsius to Kelvin."""
    return temp_c + 273.15


def _water_density_mol_per_L(temp_c: float) -> float:
    """
    Number of moles of water per liter at the given temperature.

    Implements the Tanaka et al. correlation.
    Reference: Tanaka M., Girard G., Davis R., Peuto A., Bignell N. (2001).
    """
    a1 = -3.983035
    a2 = 301.797
    a3 = 522_528.9
    a4 = 69.34881
    a5 = 999.974950

    t = temp_c
    # Density in g/L from the Tanaka correlation
    density_g_per_L = a5 * (1.0 - (t + a1) * (t + a1) * (t + a2) / a3 / (t + a4))
    # Convert to mol/L (molar mass of water = 18.0152 g/mol)
    return density_g_per_L / 18.0152


def _bimolecular(temp_c: float, G_BIMOLECULAR: float, H_BIMOLECULAR: float) -> float:
    """Bimolecular association term (kcal/mol) as a function of temperature."""
    water_density = _water_density_mol_per_L(temp_c)
    temp_k = _celcius_to_kelvin(temp_c)

    return (G_BIMOLECULAR - H_BIMOLECULAR) * temp_k / 310.15 + H_BIMOLECULAR - KB * temp_k * math.log(water_density)


def compute_assoc_energy_penalty(
    total_monomers: int,
    temp_c: float,
    G_BIMOLECULAR: float,
    H_BIMOLECULAR: float,
) -> float:
    """
    Compute the association energy penalty (kcal/mol).

    Args:
        total_monomers: Total number of monomers in the complex (>= 1).
        temp_c: Temperature in degrees Celsius.
        G_BIMOLECULAR: Empirical constant G (kcal/mol).
        H_BIMOLECULAR: Empirical constant H (kcal/mol).
        ( Nupack: G_BIMOLECULAR = 1.96, H_BIMOLECULAR = 0.20 )

    Returns:
        Association energy penalty in kcal/mol.
    """
    return _bimolecular(temp_c, G_BIMOLECULAR, H_BIMOLECULAR) * (total_monomers - 1)


class Polymer:
    """Represents a polymer as a multiset of monomers."""

    def __init__(self, monomer_counts: np.ndarray, monomers: List[Monomer], tbn: Optional["TBN"] = None):
        """
        Initialize a polymer.

        Args:
            monomer_counts: Array of monomer counts (length = number of monomers)
            monomers: List of Monomer objects in the TBN
            tbn: Optional TBN model reference for computing free energy
        """
        self.monomer_counts = monomer_counts
        self.monomers = monomers
        self.tbn = tbn
        self._concentration = None

    def get_monomers_with_counts(self) -> List[Tuple[int, Monomer]]:
        """
        Get list of (count, monomer) pairs for non-zero counts.

        Returns:
            List of tuples (count, Monomer)
        """
        result = []
        for count, monomer in zip(self.monomer_counts, self.monomers):
            if count > 0:
                result.append((int(count), monomer))
        return result

    def compute_free_energy(self, deltaG: Optional[List[float]] = None, temperature: float = 37.0) -> float:
        """
        Compute the free energy of this polymer.

        Free energy = association_energy_penalty

        Note: since the number of bonds is always maximized, 
        we are ok setting the bond free energy to zero.

        Args:
            deltaG: List of [dG_assoc, dH_assoc]. If None (default),
                   no association penalty is applied. 
            temperature: Temperature in Celsius (default: 37.0)

        Returns:
            Total free energy (bond energy + association penalty)
        """
        # Check if deltaG was explicitly provided for association penalty
        use_association_penalty = deltaG is not None
        dG_assoc = dH_assoc = 0.0
        if deltaG is not None:
            if len(deltaG) != 2:
                raise ValueError("deltaG must be [dG_assoc, dH_assoc] when provided")
            dG_assoc, dH_assoc = deltaG

        if self.tbn is None:
            raise ValueError("Cannot compute free energy without TBN model reference")

        # Bond term is ignored in this model (effectively 0)

        # Compute association energy penalty only if deltaG was explicitly provided
        if use_association_penalty:
            # Compute total number of monomers in the polymer
            total_monomers = int(np.sum(self.monomer_counts))
            # Compute association energy penalty
            assoc_penalty = compute_assoc_energy_penalty(total_monomers, temperature, dG_assoc, dH_assoc)
        else:
            # No association penalty when using default deltaG
            assoc_penalty = 0.0

        # Total free energy = association penalty only
        return assoc_penalty

    def __eq__(self, other):
        if not isinstance(other, Polymer):
            return False
        return np.array_equal(self.monomer_counts, other.monomer_counts)

    def __hash__(self):
        return hash(tuple(self.monomer_counts))


class PolymerBasisComputer:
    """Computes the polymer basis (Hilbert basis) for a TBN."""

    def __init__(
        self,
        tbn: TBN,
        normaliz_runner: Optional[NormalizRunner] = None,
        store_solver_inputs: bool = False,
        input_base_name: Optional[str] = None,
    ):
        """
        Initialize the polymer basis computer.

        Args:
            tbn: The TBN model
            normaliz_runner: Optional NormalizRunner instance (creates default if None)
            store_solver_inputs: If True, store input files for debugging
            input_base_name: Base name for stored input files
        """
        self.tbn = tbn
        self.normaliz_runner = normaliz_runner or NormalizRunner()
        self.store_solver_inputs = store_solver_inputs
        self.input_base_name = input_base_name

    def compute_polymer_basis(self) -> List[Polymer]:
        """
        Compute the polymer basis for the TBN.

        The polymer basis consists of the "unsplittable" polymers that cannot
        be decomposed into two without losing some bonding.

        Returns:
            List of Polymer objects representing the polymer basis

        Raises:
            RuntimeError: If computation fails
        """
        # Get augmented matrix A' with singleton monomers
        A_prime, n_original = self.tbn.get_augmented_matrix_for_polymer_basis()

        # Compute Hilbert basis using Normaliz
        # We want solutions to A' * x = 0 with x >= 0
        hilbert_basis_vectors = self.normaliz_runner.compute_hilbert_basis(
            A_prime,
            store_inputs=self.store_solver_inputs,
            input_base_name=self.input_base_name,
            context="polymer-basis",
        )

        if not hilbert_basis_vectors:
            raise RuntimeError("No Hilbert basis vectors found")

        # Convert Hilbert basis vectors to polymers
        # Remove entries corresponding to fake singleton monomers and remove duplicates
        polymers = []
        seen = set()

        for vector in hilbert_basis_vectors:
            # Take only the first n_original components (remove fake monomers)
            polymer_vector = vector[:n_original]

            # Create polymer object with TBN reference
            polymer = Polymer(polymer_vector, self.tbn.monomers, self.tbn)

            # Check for duplicates
            polymer_hash = hash(polymer)
            if polymer_hash not in seen:
                seen.add(polymer_hash)
                polymers.append(polymer)

        return polymers

    def save_polymer_basis(self, polymers: List[Polymer], output_file: str):
        """
        Save polymer basis to a .tbnpolys file in user-friendly format.

        Uses the .tbnpolys format specification:
        - Comments designated by "#"
        - Empty lines between polymers
        - Each polymer represented by its monomers, one per line
        - Optional multiplicity prefix "n | " where n is the count
        - Monomer shown as name (if available) or binding sites

        Args:
            polymers: List of Polymer objects
            output_file: Path to output file
        """
        from pathlib import Path

        # Convert polymers to vectors
        polymer_vectors = [polymer.monomer_counts for polymer in polymers]

        # Create writer and save
        writer = TbnpolysWriter(self.tbn)
        writer.write_polymers(
            polymer_vectors, Path(output_file), header_comment=f"Polymer basis - {len(polymers)} polymers"
        )

        print(f"Saved polymer basis with {len(polymers)} polymers to {output_file}")

    def load_cached_polymer_basis(self, polymat_file: str) -> Optional[List[Polymer]]:
        """
        Load cached polymer basis from .tbnpolymat file if matrix hash matches.

        Args:
            polymat_file: Path to .tbnpolymat file

        Returns:
            List of Polymer objects if hash matches, None otherwise
        """
        try:
            if not os.path.exists(polymat_file):
                return None

            # Check if matrix hash matches
            current_hash = self.tbn.compute_matrix_hash()
            if not check_matrix_hash(polymat_file, current_hash):
                return None

            # For backward compatibility with the old implementation:
            # Try to parse file manually to handle edge cases
            polymers = []
            has_parse_error = False

            with open(polymat_file) as f:
                for line in f:
                    line = line.strip()
                    # Skip comments and keyword lines
                    if line.startswith("#") or line.startswith("\\") or not line:
                        continue

                    parts = line.split()
                    if not parts:
                        continue

                    # Check if we have the right number of values
                    n_monomers = len(self.tbn.monomers)
                    if len(parts) < n_monomers:
                        # Wrong number of columns - skip this line
                        continue

                    # Try to parse monomer counts
                    try:
                        monomer_counts = np.array([int(x) for x in parts[:n_monomers]])
                        polymer = Polymer(monomer_counts, self.tbn.monomers, self.tbn)
                        polymers.append(polymer)
                    except ValueError:
                        # Non-numeric data - this is a parse error
                        has_parse_error = True
                        break

            # If we had parse errors (non-numeric data), return None
            if has_parse_error:
                return None

            # Otherwise return the polymers (could be empty list)
            return polymers

        except Exception:
            # If any error occurs during loading, return None to recompute
            return None

    def save_tbnpolymat(
        self,
        polymers: List[Polymer],
        output_file: str,
        compute_free_energies: bool = True,
        compute_concentrations: bool = True,
        concentration_runner: Optional[object] = None,
        verbose: bool = False,
        parameters: Optional[dict] = None,
        deltaG: Optional[List[float]] = None,
        temperature: float = 37.0,
    ):
        """
        Save polymer basis to .tbnpolymat file format.

        Format: Each row has monomer counts, optionally free energy, optionally concentration.
        Polymers are sorted by decreasing concentration (if available).

        Args:
            polymers: List of Polymer objects
            output_file: Path to output .tbnpolymat file
            compute_free_energies: Whether to compute and include free energies
            compute_concentrations: Whether to compute and include concentrations
            concentration_runner: Optional COFFEERunner or NupackRunner instance for concentration computation
            verbose: Whether to enable verbose output
            parameters: Optional dictionary of parameters used for parametrized .tbn files
            deltaG: Optional [dG_assoc, dH_assoc] association parameters
            temperature: Temperature in Celsius (default: 37.0)
        """
        # Determine what to compute
        has_monomer_concentrations = self.tbn.concentrations is not None
        include_free_energies = compute_free_energies
        include_concentrations = (
            compute_concentrations
            and compute_free_energies  # Can't compute concentrations without free energies
            and has_monomer_concentrations
        )

        # Compute concentrations if requested and possible
        polymer_concentrations = None
        if include_concentrations:
            if concentration_runner is None:
                concentration_runner = COFFEERunner()
            try:
                polymer_concentrations = concentration_runner.compute_equilibrium_concentrations(
                    polymers, self.tbn, deltaG=deltaG, temperature=temperature
                )
                if verbose:
                    print("Equilibrium concentrations computed")
                # Attach concentrations to polymers for sorting
                for polymer, conc in zip(polymers, polymer_concentrations):
                    polymer._concentration = conc
            except Exception as e:
                print(f"Warning: Could not compute concentrations: {e}")
                include_concentrations = False

        # Sort polymers by concentration if available
        if include_concentrations and polymer_concentrations is not None:
            # Sort by concentration in descending order
            sorted_indices = np.argsort(-polymer_concentrations)
            sorted_polymers = [polymers[i] for i in sorted_indices]
            sorted_concentrations = polymer_concentrations[sorted_indices]
        else:
            sorted_polymers = polymers
            sorted_concentrations = None

        # Prepare polymer data for PolymatData
        polymer_arrays = [polymer.monomer_counts for polymer in sorted_polymers]

        # Compute free energies if requested
        free_energies = None
        if include_free_energies:
            free_energies = np.array([polymer.compute_free_energy(deltaG, temperature) for polymer in sorted_polymers])

        # Create PolymatData object
        polymat_data = PolymatData(
            polymers=polymer_arrays,
            n_monomers=len(self.tbn.monomers),
            n_polymers=len(sorted_polymers),
            matrix_hash=self.tbn.compute_matrix_hash(),
            free_energies=free_energies,
            concentrations=sorted_concentrations,
            concentration_units=get_unit_display_name(self.tbn.concentration_units) if include_concentrations else None,
            has_free_energies=include_free_energies,
            has_concentrations=include_concentrations,
            parameters=parameters,
        )

        # Write using PolymatWriter
        writer = PolymatWriter(output_file)
        writer.write(polymat_data, from_molar, self.tbn.concentration_units)
