import numpy as np
import os
from typing import List, Tuple, Optional
from .model import TBN, Monomer
from .normaliz import NormalizRunner
from .coffee import COFFEERunner
from .units import from_molar, get_unit_display_name
from .polymat_io import PolymatData, PolymatReader, PolymatWriter, check_matrix_hash


class Polymer:
    """Represents a polymer as a multiset of monomers."""
    
    def __init__(self, monomer_counts: np.ndarray, monomers: List[Monomer], tbn: Optional['TBN'] = None):
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
        self._free_energy = None
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
    
    def compute_free_energy(self) -> float:
        """
        Compute the free energy of this polymer.
        
        Free energy = -[number of bonds in polymer]
        Number of bonds = (Sum[|A|.x] - Sum[A.x])/2
        
        Where:
        - |A| is the matrix A with absolute values
        - x is the polymer's monomer count vector
        - Sum[v] sums all elements of vector v
        - We divide by 2 because each bond involves exactly 2 binding sites
        
        Returns:
            Free energy (negative number of bonds)
        """
        if self._free_energy is not None:
            return self._free_energy
        
        if self.tbn is None:
            raise ValueError("Cannot compute free energy without TBN model reference")
        
        # Get matrix A
        A = self.tbn.matrix_A
        
        # Compute |A| * x (total binding sites, excluding self-binding within monomers)
        abs_A = np.abs(A)
        total_binding_sites = np.sum(abs_A @ self.monomer_counts)
        
        # Compute A * x (excess of unstar binding sites)
        excess_unstar = np.sum(A @ self.monomer_counts)
        
        # Number of bonds = (total binding sites - excess unstar) / 2
        # Divide by 2 because each bond involves exactly 2 binding sites
        num_bonds = (total_binding_sites - excess_unstar) / 2
        
        # Free energy = -number of bonds
        self._free_energy = -num_bonds
        
        return self._free_energy
    
    def __eq__(self, other):
        if not isinstance(other, Polymer):
            return False
        return np.array_equal(self.monomer_counts, other.monomer_counts)
    
    def __hash__(self):
        return hash(tuple(self.monomer_counts))


class PolymerBasisComputer:
    """Computes the polymer basis (Hilbert basis) for a TBN."""
    
    def __init__(self, tbn: TBN, normaliz_runner: Optional[NormalizRunner] = None):
        """
        Initialize the polymer basis computer.
        
        Args:
            tbn: The TBN model
            normaliz_runner: Optional NormalizRunner instance (creates default if None)
        """
        self.tbn = tbn
        self.normaliz_runner = normaliz_runner or NormalizRunner()
    
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
        hilbert_basis_vectors = self.normaliz_runner.compute_hilbert_basis(A_prime)
        
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
        Save polymer basis to a text file in user-friendly format.
        
        Format:
        - Empty line between polymers
        - Each polymer represented by its monomers, one per line
        - Line starts with "n | " where n is the multiplicity
        - Monomer shown as name (if available) or binding sites
        
        Args:
            polymers: List of Polymer objects
            output_file: Path to output file
        """
        with open(output_file, 'w') as f:
            f.write(f"# Polymer basis - {len(polymers)} polymers\n")
            f.write("#\n")
            
            for i, polymer in enumerate(polymers):
                if i > 0:
                    f.write("\n")  # Empty line between polymers
                
                f.write(f"# Polymer {i+1}\n")
                
                # Get monomers with their counts
                monomers_with_counts = polymer.get_monomers_with_counts()
                
                for count, monomer in monomers_with_counts:
                    # Format: "count | monomer_representation"
                    if monomer.name:
                        monomer_str = monomer.name
                    else:
                        monomer_str = monomer.get_binding_sites_str()
                    
                    f.write(f"{count} | {monomer_str}\n")
        
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
            
            with open(polymat_file, 'r') as f:
                for line in f:
                    line = line.strip()
                    # Skip comments and keyword lines
                    if line.startswith('#') or line.startswith('\\') or not line:
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
        coffee_runner: Optional[COFFEERunner] = None,
        verbose: bool = False
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
            coffee_runner: Optional COFFEERunner instance for concentration computation
        """
        # Determine what to compute
        has_monomer_concentrations = self.tbn.concentrations is not None
        include_free_energies = compute_free_energies
        include_concentrations = (
            compute_concentrations and 
            compute_free_energies and  # Can't compute concentrations without free energies
            has_monomer_concentrations
        )
        
        # Compute concentrations if requested and possible
        polymer_concentrations = None
        if include_concentrations:
            if coffee_runner is None:
                coffee_runner = COFFEERunner()
            try:
                polymer_concentrations = coffee_runner.compute_equilibrium_concentrations(
                    polymers, self.tbn
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
            free_energies = np.array([polymer.compute_free_energy() for polymer in sorted_polymers])
        
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
            has_concentrations=include_concentrations
        )
        
        # Write using PolymatWriter
        writer = PolymatWriter(output_file)
        writer.write(polymat_data, from_molar, self.tbn.concentration_units)