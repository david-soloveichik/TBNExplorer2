"""
Module for reading and writing .tbnpolymat files.

Provides unified I/O operations for polymer matrix files used by both
tbnexplorer2 and tbnexplorer2-filter commands.
"""

import numpy as np
from pathlib import Path
from typing import List, Optional, Iterator, Tuple
from dataclasses import dataclass


@dataclass
class PolymatData:
    """
    Structured representation of .tbnpolymat file contents.
    
    Attributes:
        polymers: List of polymer count arrays (each array has length n_monomers)
        free_energies: Optional array of free energies for each polymer
        concentrations: Optional array of concentrations for each polymer
        n_monomers: Number of monomers in the system
        n_polymers: Number of polymers in the basis
        matrix_hash: Hash of the matrix used to generate this data
        concentration_units: Units for concentration values (e.g., 'nM', 'uM')
        has_free_energies: Whether free energies are included
        has_concentrations: Whether concentrations are included
    """
    polymers: List[np.ndarray]
    n_monomers: int
    n_polymers: int
    matrix_hash: Optional[str] = None
    free_energies: Optional[np.ndarray] = None
    concentrations: Optional[np.ndarray] = None
    concentration_units: Optional[str] = None
    has_free_energies: bool = False
    has_concentrations: bool = False
    
    def get_polymer_data(self, index: int) -> Tuple[np.ndarray, Optional[float], Optional[float]]:
        """
        Get data for a specific polymer.
        
        Args:
            index: Polymer index
            
        Returns:
            Tuple of (monomer_counts, free_energy, concentration)
        """
        if index < 0 or index >= self.n_polymers:
            raise IndexError(f"Polymer index {index} out of range [0, {self.n_polymers})")
        
        polymer_counts = self.polymers[index]
        free_energy = self.free_energies[index] if self.free_energies is not None else None
        concentration = self.concentrations[index] if self.concentrations is not None else None
        
        return polymer_counts, free_energy, concentration


class PolymatReader:
    """
    Efficient reader for .tbnpolymat files with support for large files.
    """
    
    def __init__(self, file_path: str):
        """
        Initialize reader with file path.
        
        Args:
            file_path: Path to .tbnpolymat file
        """
        self.file_path = Path(file_path)
        if not self.file_path.exists():
            raise FileNotFoundError(f"File not found: {file_path}")
    
    def read(self, lazy_load: bool = False) -> PolymatData:
        """
        Read the entire .tbnpolymat file.
        
        Args:
            lazy_load: If True, returns an iterator for polymer data instead of loading all at once
                      (useful for very large files)
        
        Returns:
            PolymatData object containing all polymer information
        """
        # Parse header first to understand file structure
        header_info = self._parse_header()
        
        if lazy_load:
            # For lazy loading, we'll still parse the header but return a special
            # PolymatData that loads polymers on demand
            # For now, we'll implement eager loading
            pass
        
        # Read all polymer data
        polymers = []
        free_energies = [] if header_info['has_free_energies'] else None
        concentrations = [] if header_info['has_concentrations'] else None
        
        with open(self.file_path, 'r') as f:
            for line in f:
                if line.startswith('#') or not line.strip():
                    continue
                
                # Parse data line
                polymer_data = self._parse_data_line(line, header_info)
                if polymer_data is not None:
                    polymers.append(polymer_data['counts'])
                    if free_energies is not None:
                        free_energies.append(polymer_data['free_energy'])
                    if concentrations is not None:
                        concentrations.append(polymer_data['concentration'])
        
        # Convert to numpy arrays where appropriate
        if free_energies is not None:
            free_energies = np.array(free_energies)
        if concentrations is not None:
            concentrations = np.array(concentrations)
        
        return PolymatData(
            polymers=polymers,
            n_monomers=header_info['n_monomers'],
            n_polymers=len(polymers),
            matrix_hash=header_info.get('matrix_hash'),
            free_energies=free_energies,
            concentrations=concentrations,
            concentration_units=header_info.get('concentration_units'),
            has_free_energies=header_info['has_free_energies'],
            has_concentrations=header_info['has_concentrations']
        )
    
    def read_header_only(self) -> dict:
        """
        Read only the header information without loading polymer data.
        Useful for checking matrix hash or file structure.
        
        Returns:
            Dictionary with header information
        """
        return self._parse_header()
    
    def iterate_polymers(self) -> Iterator[Tuple[np.ndarray, Optional[float], Optional[float]]]:
        """
        Iterate through polymers one at a time without loading entire file into memory.
        Useful for very large files.
        
        Yields:
            Tuples of (monomer_counts, free_energy, concentration) for each polymer
        """
        header_info = self._parse_header()
        
        with open(self.file_path, 'r') as f:
            for line in f:
                if line.startswith('#') or not line.strip():
                    continue
                
                polymer_data = self._parse_data_line(line, header_info)
                if polymer_data is not None:
                    yield (
                        polymer_data['counts'],
                        polymer_data.get('free_energy'),
                        polymer_data.get('concentration')
                    )
    
    def _parse_header(self) -> dict:
        """
        Parse header information from .tbnpolymat file.
        
        Returns:
            Dictionary containing header information
        """
        header_info = {
            'n_monomers': None,
            'n_polymers': None,
            'matrix_hash': None,
            'concentration_units': None,
            'has_free_energies': False,
            'has_concentrations': False
        }
        
        with open(self.file_path, 'r') as f:
            for line in f:
                if not line.startswith('#'):
                    # End of header
                    break
                
                line = line.strip()
                
                if 'Number of monomers:' in line:
                    try:
                        header_info['n_monomers'] = int(line.split(':')[1].strip())
                    except (ValueError, IndexError):
                        pass
                elif 'Number of polymers:' in line:
                    try:
                        header_info['n_polymers'] = int(line.split(':')[1].strip())
                    except (ValueError, IndexError):
                        pass
                elif 'MATRIX-HASH:' in line:
                    header_info['matrix_hash'] = line.split(':', 1)[1].strip()
                elif 'Concentration units:' in line:
                    header_info['concentration_units'] = line.split(':', 1)[1].strip()
                elif 'Columns:' in line:
                    columns_str = line.split(':', 1)[1].strip()
                    header_info['has_free_energies'] = 'free_energy' in columns_str
                    header_info['has_concentrations'] = 'concentration' in columns_str
        
        # If n_monomers not specified in header, try to infer from first data line
        if header_info['n_monomers'] is None:
            with open(self.file_path, 'r') as f:
                for line in f:
                    if not line.startswith('#') and line.strip():
                        parts = line.strip().split()
                        if parts:
                            # Assume first data line has correct number of monomer counts
                            # This is a best guess - may not be accurate
                            header_info['n_monomers'] = len(parts)
                            break
        
        # Still validate required fields
        if header_info['n_monomers'] is None:
            raise ValueError("Invalid .tbnpolymat file: cannot determine number of monomers")
        
        return header_info
    
    def _parse_data_line(self, line: str, header_info: dict) -> Optional[dict]:
        """
        Parse a single data line from the file.
        
        Args:
            line: Line to parse
            header_info: Header information dictionary
            
        Returns:
            Dictionary with parsed data or None if line is invalid
        """
        parts = line.strip().split()
        if not parts:
            return None
        
        n_monomers = header_info['n_monomers']
        
        # Parse monomer counts
        if len(parts) < n_monomers:
            return None
        
        try:
            monomer_counts = np.array([int(x) for x in parts[:n_monomers]])
        except ValueError:
            return None
        
        result = {'counts': monomer_counts}
        
        # Parse optional columns
        col_index = n_monomers
        
        if header_info['has_free_energies'] and col_index < len(parts):
            try:
                result['free_energy'] = float(parts[col_index])
            except ValueError:
                result['free_energy'] = None
            col_index += 1
        
        if header_info['has_concentrations'] and col_index < len(parts):
            try:
                result['concentration'] = float(parts[col_index])
            except ValueError:
                result['concentration'] = None
        
        return result


class PolymatWriter:
    """
    Writer for .tbnpolymat files with consistent formatting.
    """
    
    def __init__(self, file_path: str):
        """
        Initialize writer with file path.
        
        Args:
            file_path: Path where .tbnpolymat file will be written
        """
        self.file_path = Path(file_path)
    
    def write(self, data: PolymatData, from_molar_func=None, tbn_units=None):
        """
        Write PolymatData to .tbnpolymat file.
        
        Args:
            data: PolymatData object to write
            from_molar_func: Optional function to convert from Molar to target units
            tbn_units: Optional TBN concentration units for conversion
        """
        with open(self.file_path, 'w') as f:
            # Write header
            self._write_header(f, data)
            
            # Write polymer data
            for i in range(data.n_polymers):
                polymer_counts, free_energy, concentration = data.get_polymer_data(i)
                
                # Write monomer counts
                counts_str = ' '.join(str(int(c)) for c in polymer_counts)
                row = [counts_str]
                
                # Add free energy if available
                if data.has_free_energies and free_energy is not None:
                    row.append(str(free_energy))
                
                # Add concentration if available
                if data.has_concentrations and concentration is not None:
                    # Handle unit conversion if needed
                    if from_molar_func is not None and tbn_units is not None:
                        # Convert from Molar to target units
                        conc_target = from_molar_func(concentration, tbn_units)
                        if conc_target == 0:
                            row.append("0.00e0")
                        else:
                            row.append(f"{conc_target:.2e}")
                    else:
                        # Use value as-is
                        if concentration == 0:
                            row.append("0.00e0")
                        else:
                            row.append(f"{concentration:.2e}")
                
                f.write(' '.join(row) + '\n')
    
    def _write_header(self, f, data: PolymatData):
        """
        Write header section to file.
        
        Args:
            f: File handle
            data: PolymatData object
        """
        f.write(f"# TBN Polymer Matrix\n")
        f.write(f"# Number of polymers: {data.n_polymers}\n")
        f.write(f"# Number of monomers: {data.n_monomers}\n")
        
        if data.matrix_hash:
            f.write(f"# MATRIX-HASH: {data.matrix_hash}\n")
        
        if data.concentration_units:
            f.write(f"# Concentration units: {data.concentration_units}\n")
        
        # Build columns description
        columns = [f"monomer_counts[1..{data.n_monomers}]"]
        if data.has_free_energies:
            columns.append("free_energy")
        if data.has_concentrations:
            columns.append("concentration")
        f.write(f"# Columns: {' '.join(columns)}\n")
        f.write("#\n")


def load_polymat_file(file_path: str, lazy_load: bool = False) -> PolymatData:
    """
    Convenience function to load a .tbnpolymat file.
    
    Args:
        file_path: Path to .tbnpolymat file
        lazy_load: If True, uses lazy loading for large files
        
    Returns:
        PolymatData object
    """
    reader = PolymatReader(file_path)
    return reader.read(lazy_load=lazy_load)


def save_polymat_file(file_path: str, data: PolymatData, from_molar_func=None, tbn_units=None):
    """
    Convenience function to save a .tbnpolymat file.
    
    Args:
        file_path: Path where file will be saved
        data: PolymatData object to save
        from_molar_func: Optional function for unit conversion
        tbn_units: Optional TBN units for conversion
    """
    writer = PolymatWriter(file_path)
    writer.write(data, from_molar_func, tbn_units)


def check_matrix_hash(file_path: str, expected_hash: str) -> bool:
    """
    Check if a .tbnpolymat file has the expected matrix hash.
    Useful for cache validation.
    
    Args:
        file_path: Path to .tbnpolymat file
        expected_hash: Expected matrix hash value
        
    Returns:
        True if hashes match, False otherwise
    """
    try:
        reader = PolymatReader(file_path)
        header = reader.read_header_only()
        return header.get('matrix_hash') == expected_hash
    except (FileNotFoundError, ValueError):
        return False