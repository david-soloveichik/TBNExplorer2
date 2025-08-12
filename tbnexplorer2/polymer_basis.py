import numpy as np
from typing import List, Tuple, Optional
from .model import TBN, Monomer
from .normaliz import NormalizRunner
from .coffee import COFFEERunner


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
    
    def save_tbnpolymat(
        self, 
        polymers: List[Polymer], 
        output_file: str,
        compute_free_energies: bool = True,
        compute_concentrations: bool = True,
        coffee_runner: Optional[COFFEERunner] = None
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
            if polymer_concentrations is not None:
                sorted_concentrations = polymer_concentrations[sorted_indices]
        else:
            sorted_polymers = polymers
            sorted_concentrations = None
        
        # Write the file
        with open(output_file, 'w') as f:
            # Write header comments
            f.write(f"# TBN Polymer Matrix\n")
            f.write(f"# Number of polymers: {len(sorted_polymers)}\n")
            f.write(f"# Number of monomers: {len(self.tbn.monomers)}\n")
            
            # Indicate what's included
            computation_status = []
            if not has_monomer_concentrations:
                computation_status.append("no monomer concentrations provided")
            if not include_free_energies:
                computation_status.append("free energies not computed")
            if has_monomer_concentrations and not include_concentrations:
                computation_status.append("equilibrium concentrations not computed")
            
            if computation_status:
                f.write(f"# Partial computation: {', '.join(computation_status)}\n")
            
            # Write column description
            columns = ["monomer_counts[1..{}]".format(len(self.tbn.monomers))]
            if include_free_energies:
                columns.append("free_energy")
            if include_concentrations:
                columns.append("concentration")
            f.write(f"# Columns: {' '.join(columns)}\n")
            f.write("#\n")
            
            # Write polymer data
            for i, polymer in enumerate(sorted_polymers):
                # Write monomer counts
                counts_str = ' '.join(str(int(c)) for c in polymer.monomer_counts)
                row = [counts_str]
                
                # Add free energy if requested
                if include_free_energies:
                    free_energy = polymer.compute_free_energy()
                    row.append(str(free_energy))
                
                # Add concentration if available
                if include_concentrations and sorted_concentrations is not None:
                    conc = sorted_concentrations[i]
                    # Format scientific notation properly
                    if conc == 0:
                        row.append("0.00e0")
                    else:
                        row.append(f"{conc:.2e}")
                
                f.write(' '.join(row) + '\n')