import numpy as np
import hashlib
from typing import List, Optional, Dict, Tuple
from dataclasses import dataclass
from .units import to_molar


@dataclass
class BindingSite:
    """Represents a binding site in the TBN model."""
    name: str
    is_star: bool
    
    def __str__(self):
        return f"{self.name}*" if self.is_star else self.name
    
    def __repr__(self):
        return f"BindingSite('{self.name}', star={self.is_star})"
    
    def __eq__(self, other):
        if not isinstance(other, BindingSite):
            return False
        return self.name == other.name and self.is_star == other.is_star
    
    def __hash__(self):
        return hash((self.name, self.is_star))


@dataclass
class Monomer:
    """Represents a monomer in the TBN model."""
    name: Optional[str]
    binding_sites: List[BindingSite]
    concentration: Optional[float]
    original_line: str
    
    def to_vector(self, binding_site_index: Dict[str, int]) -> np.ndarray:
        """
        Convert monomer to integer vector representation.
        
        Star binding sites contribute -1, unstar binding sites contribute +1.
        
        Args:
            binding_site_index: Dictionary mapping binding site names to indices
            
        Returns:
            numpy array representing the monomer
        """
        vector = np.zeros(len(binding_site_index), dtype=int)
        
        for site in self.binding_sites:
            idx = binding_site_index[site.name]
            if site.is_star:
                vector[idx] -= 1
            else:
                vector[idx] += 1
        
        return vector
    
    def get_binding_sites_str(self) -> str:
        """Get the binding sites as a string in original order."""
        return ' '.join(str(site) for site in self.binding_sites)
    
    def __str__(self):
        sites_str = self.get_binding_sites_str()
        if self.name:
            return f"{self.name}: {sites_str}"
        return sites_str


class TBN:
    """Represents a complete TBN (Thermodynamics of Binding Networks) model."""
    
    def __init__(self, monomers: List[Monomer], binding_site_index: Dict[str, int], concentration_units: Optional[str] = None):
        """
        Initialize a TBN model.
        
        Args:
            monomers: List of Monomer objects
            binding_site_index: Dictionary mapping binding site names to indices
            concentration_units: Units for input/output concentrations or None if no concentrations
        """
        self.monomers = monomers
        self.binding_site_index = binding_site_index
        self.concentration_units = concentration_units
        self._matrix_A = None
        self._concentrations = None
        self._concentrations_molar = None
    
    @property
    def matrix_A(self) -> np.ndarray:
        """
        Get the matrix A where columns are monomer vectors.
        
        Returns:
            numpy array of shape (n_binding_sites, n_monomers)
        """
        if self._matrix_A is None:
            if not self.monomers:
                # Handle empty monomers case
                self._matrix_A = np.zeros((len(self.binding_site_index), 0), dtype=int)
            else:
                vectors = [m.to_vector(self.binding_site_index) for m in self.monomers]
                self._matrix_A = np.column_stack(vectors)
        return self._matrix_A
    
    @property
    def concentrations(self) -> Optional[np.ndarray]:
        """
        Get the concentration vector in Molar units if all monomers have concentrations.
        
        Returns:
            numpy array of concentrations in Molar or None
        """
        if self._concentrations_molar is None:
            if all(m.concentration is not None for m in self.monomers) and self.concentration_units is not None:
                # Convert from input units to Molar for internal processing
                original_concentrations = np.array([m.concentration for m in self.monomers])
                self._concentrations_molar = to_molar(original_concentrations, self.concentration_units)
        return self._concentrations_molar
    
    @property
    def concentrations_original_units(self) -> Optional[np.ndarray]:
        """
        Get the concentration vector in original input units if all monomers have concentrations.
        
        Returns:
            numpy array of concentrations in original units or None
        """
        if self._concentrations is None:
            if all(m.concentration is not None for m in self.monomers):
                self._concentrations = np.array([m.concentration for m in self.monomers])
        return self._concentrations
    
    def check_star_limiting(self) -> Tuple[bool, Optional[str]]:
        """
        Check if the TBN satisfies the star-limiting restriction.
        
        The TBN is star-limited if for every binding site there is at least 
        as much unstar as star (totalled over all monomers).
        
        Returns:
            Tuple of (is_valid, error_message or None)
        """
        if self.concentrations is not None:
            # Use provided concentrations
            c = self.concentrations
        else:
            # Use unit concentrations
            c = np.ones(len(self.monomers))
        
        # Compute A * c
        binding_site_excesses = self.matrix_A @ c
        
        # Check if all components are non-negative
        if np.all(binding_site_excesses >= 0):
            return True, None
        
        # Find problematic binding sites
        problematic_sites = []
        for i, excess in enumerate(binding_site_excesses):
            if excess < 0:
                site_name = None
                for name, idx in self.binding_site_index.items():
                    if idx == i:
                        site_name = name
                        break
                problematic_sites.append((site_name, excess))
        
        error_msg = "TBN is not star-limited. Binding sites with negative excess:\n"
        for site, excess in problematic_sites:
            error_msg += f"  {site}: {excess:.2f}\n"
        
        return False, error_msg
    
    def get_augmented_matrix_for_polymer_basis(self) -> Tuple[np.ndarray, int]:
        """
        Get the augmented matrix A' for polymer basis computation.
        
        For every binding site x that doesn't have a singleton {x*} monomer,
        add a column with all zeros and -1 in the position for x.
        
        Returns:
            Tuple of (augmented matrix A', number of original monomers)
        """
        A = self.matrix_A.copy()
        n_original = A.shape[1]
        
        # Check which binding sites need singleton monomers
        additional_columns = []
        
        for site_name, site_idx in self.binding_site_index.items():
            # Check if there's already a singleton {x*} monomer
            need_singleton = True
            
            for col_idx in range(n_original):
                col = A[:, col_idx]
                # Check if this column is all zeros except -1 at site_idx
                if col[site_idx] == -1 and np.sum(np.abs(col)) == 1:
                    need_singleton = False
                    break
            
            if need_singleton:
                # Create singleton column
                singleton_col = np.zeros(len(self.binding_site_index), dtype=int)
                singleton_col[site_idx] = -1
                additional_columns.append(singleton_col)
        
        # Augment matrix if needed
        if additional_columns:
            A_prime = np.column_stack([A] + additional_columns)
        else:
            A_prime = A
        
        return A_prime, n_original
    
    def compute_matrix_hash(self) -> str:
        """
        Compute SHA256 hash of the matrix A for caching purposes.
        
        The hash is computed from the matrix A in a deterministic way to enable
        caching of polymer basis computations.
        
        Returns:
            Hexadecimal string representation of the hash
        """
        # Convert matrix to bytes in a deterministic way
        # Use the matrix shape and flattened contents
        matrix_data = self.matrix_A.tobytes()
        shape_data = np.array(self.matrix_A.shape).tobytes()
        
        # Combine shape and matrix data
        combined_data = shape_data + matrix_data
        
        # Compute SHA256 hash
        hash_obj = hashlib.sha256(combined_data)
        return hash_obj.hexdigest()
    
    def __str__(self):
        n_sites = len(self.binding_site_index)
        n_monomers = len(self.monomers)
        has_conc = self.concentrations is not None
        
        return (f"TBN Model:\n"
                f"  Binding sites: {n_sites}\n"
                f"  Monomers: {n_monomers}\n"
                f"  Has concentrations: {has_conc}")
    
    def __repr__(self):
        return f"TBN(monomers={len(self.monomers)}, binding_sites={len(self.binding_site_index)})"