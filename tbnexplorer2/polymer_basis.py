import numpy as np
from typing import List, Tuple, Optional
from .model import TBN, Monomer
from .normaliz import NormalizRunner


class Polymer:
    """Represents a polymer as a multiset of monomers."""
    
    def __init__(self, monomer_counts: np.ndarray, monomers: List[Monomer]):
        """
        Initialize a polymer.
        
        Args:
            monomer_counts: Array of monomer counts (length = number of monomers)
            monomers: List of Monomer objects in the TBN
        """
        self.monomer_counts = monomer_counts
        self.monomers = monomers
    
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
            
            # Create polymer object
            polymer = Polymer(polymer_vector, self.tbn.monomers)
            
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