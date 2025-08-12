"""
Filter module for TBN Explorer 2.

Filters polymers from .tbnpolymat files based on monomer name criteria.
"""

import numpy as np
from pathlib import Path
from typing import List, Dict, Tuple, Optional
from .parser import TBNParser
from .model import TBN, Monomer


class PolymerFilter:
    """Filters polymers based on monomer name criteria."""
    
    def __init__(self, tbn_file: str):
        """
        Initialize the filter with a TBN file.
        
        Args:
            tbn_file: Path to the .tbn file
        """
        self.tbn_file = Path(tbn_file)
        self.polymat_file = self._infer_polymat_file()
        
        # Parse TBN file to get monomer information
        self.monomers, self.binding_site_index = TBNParser.parse_file(str(self.tbn_file))
        
        # Load polymer data from .tbnpolymat file
        self.polymer_data = self._load_polymat_file()
        
    def _infer_polymat_file(self) -> Path:
        """
        Infer the .tbnpolymat file path from the .tbn file.
        
        Returns:
            Path to the .tbnpolymat file
            
        Raises:
            FileNotFoundError: If the .tbnpolymat file doesn't exist
        """
        base_name = self.tbn_file.stem
        polymat_file = self.tbn_file.parent / f"{base_name}.tbnpolymat"
        
        if not polymat_file.exists():
            raise FileNotFoundError(f"Cannot find polymer matrix file: {polymat_file}")
        
        return polymat_file
    
    def _load_polymat_file(self) -> Dict:
        """
        Load and parse the .tbnpolymat file.
        
        Returns:
            Dictionary containing polymer data with keys:
            - 'polymers': List of polymer count arrays
            - 'free_energies': Optional array of free energies
            - 'concentrations': Optional array of concentrations
            - 'concentration_units': Optional concentration units string
            - 'has_free_energies': Boolean indicating if free energies are present
            - 'has_concentrations': Boolean indicating if concentrations are present
        """
        data = {
            'polymers': [],
            'free_energies': None,
            'concentrations': None,
            'concentration_units': None,
            'has_free_energies': False,
            'has_concentrations': False
        }
        
        free_energies = []
        concentrations = []
        
        with open(self.polymat_file, 'r') as f:
            # Parse header to determine what columns are present
            for line in f:
                if line.startswith('# Concentration units:'):
                    # Extract units like "nanoMolar (nM)" -> "nM"
                    units_str = line.split(':', 1)[1].strip()
                    if '(' in units_str and ')' in units_str:
                        data['concentration_units'] = units_str.split('(')[1].rstrip(')')
                    else:
                        # Fallback to parsing the first word
                        data['concentration_units'] = units_str.split()[0]
                elif line.startswith('# Columns:'):
                    columns_str = line.split(':', 1)[1].strip()
                    data['has_free_energies'] = 'free_energy' in columns_str
                    data['has_concentrations'] = 'concentration' in columns_str
                elif not line.startswith('#'):
                    # Data line
                    parts = line.strip().split()
                    if not parts:
                        continue
                    
                    # First n_monomers values are monomer counts
                    n_monomers = len(self.monomers)
                    monomer_counts = np.array([int(x) for x in parts[:n_monomers]])
                    data['polymers'].append(monomer_counts)
                    
                    # Check for additional columns
                    col_index = n_monomers
                    if data['has_free_energies'] and col_index < len(parts):
                        free_energies.append(float(parts[col_index]))
                        col_index += 1
                    
                    if data['has_concentrations'] and col_index < len(parts):
                        # Parse scientific notation properly
                        conc_str = parts[col_index]
                        concentrations.append(float(conc_str))
        
        # Convert lists to arrays
        if free_energies:
            data['free_energies'] = np.array(free_energies)
        if concentrations:
            data['concentrations'] = np.array(concentrations)
        
        return data
    
    def filter_by_monomers(self, monomer_names: List[str], percent_limit: Optional[float] = None) -> List[Tuple[int, np.ndarray, Optional[float], Optional[float]]]:
        """
        Filter polymers containing all specified monomers.
        
        Args:
            monomer_names: List of monomer names to filter by (can have duplicates for multiplicity)
            percent_limit: Optional percentage limit (0-100) for filtering by concentration
            
        Returns:
            List of tuples (polymer_index, monomer_counts, free_energy, concentration)
            sorted by decreasing concentration (if available) or by polymer index
        """
        # Count required multiplicity for each monomer name
        required_counts = {}
        for name in monomer_names:
            required_counts[name] = required_counts.get(name, 0) + 1
        
        # Create mapping from monomer name/binding sites to indices
        monomer_name_to_indices = {}
        for i, monomer in enumerate(self.monomers):
            # Use monomer name if available, otherwise use binding sites string
            identifier = monomer.name if monomer.name else monomer.get_binding_sites_str()
            if identifier not in monomer_name_to_indices:
                monomer_name_to_indices[identifier] = []
            monomer_name_to_indices[identifier].append(i)
        
        # Check if all required monomer names exist
        for name in required_counts:
            if name not in monomer_name_to_indices:
                # Return empty list if monomer name doesn't exist
                return []
        
        # Calculate total concentration if available
        total_concentration = None
        if self.polymer_data['has_concentrations'] and self.polymer_data['concentrations'] is not None:
            total_concentration = np.sum(self.polymer_data['concentrations'])
        
        # Filter polymers
        matching_polymers = []
        for i, polymer_counts in enumerate(self.polymer_data['polymers']):
            # Check if polymer contains all required monomers with correct multiplicity
            matches = True
            for monomer_name, required_count in required_counts.items():
                # Sum counts for all monomers with this name
                actual_count = sum(polymer_counts[idx] for idx in monomer_name_to_indices[monomer_name])
                if actual_count < required_count:
                    matches = False
                    break
            
            if matches:
                free_energy = None
                if self.polymer_data['free_energies'] is not None:
                    free_energy = self.polymer_data['free_energies'][i]
                
                concentration = None
                if self.polymer_data['concentrations'] is not None:
                    concentration = self.polymer_data['concentrations'][i]
                
                # Apply percent limit if specified
                if percent_limit is not None and concentration is not None and total_concentration is not None:
                    concentration_percent = (concentration / total_concentration) * 100
                    if concentration_percent < percent_limit:
                        continue
                
                matching_polymers.append((i, polymer_counts, free_energy, concentration))
        
        # Sort by concentration (descending) if available
        if self.polymer_data['has_concentrations']:
            matching_polymers.sort(key=lambda x: x[3] if x[3] is not None else 0, reverse=True)
        
        return matching_polymers
    
    def format_output(self, filtered_polymers: List[Tuple[int, np.ndarray, Optional[float], Optional[float]]], 
                     monomer_names: List[str], percent_limit: Optional[float] = None) -> str:
        """
        Format filtered polymers for user-friendly output.
        
        Args:
            filtered_polymers: List of filtered polymer tuples
            monomer_names: Original list of monomer names used for filtering
            percent_limit: Percent limit used for filtering (if any)
            
        Returns:
            Formatted string output
        """
        output_lines = []
        
        # Header
        output_lines.append(f"# Filtered polymers containing: {' '.join(monomer_names)}")
        if percent_limit is not None:
            output_lines.append(f"# Percent limit: {percent_limit}%")
        output_lines.append(f"# Number of matching polymers: {len(filtered_polymers)}")
        
        # Calculate total concentration stats if available
        if self.polymer_data['has_concentrations'] and self.polymer_data['concentrations'] is not None:
            total_concentration = np.sum(self.polymer_data['concentrations'])
            matching_concentration = sum(p[3] for p in filtered_polymers if p[3] is not None)
            percentage = (matching_concentration / total_concentration * 100) if total_concentration > 0 else 0
            output_lines.append(f"# Total concentration fraction: {percentage:.2f}%")
            if self.polymer_data['concentration_units']:
                output_lines.append(f"# Concentration units: {self.polymer_data['concentration_units']}")
        
        output_lines.append("#")
        
        # Format each polymer
        for idx, (polymer_idx, monomer_counts, free_energy, concentration) in enumerate(filtered_polymers, 1):
            output_lines.append(f"# Polymer {idx}")
            
            # Show polymer composition
            for count, monomer in zip(monomer_counts, self.monomers):
                if count > 0:
                    sites_str = monomer.get_binding_sites_str()
                    if monomer.name:
                        # Named monomer
                        output_lines.append(f"{count} | {monomer.name}")
                    else:
                        # Unnamed monomer - show binding sites
                        output_lines.append(f"{count} | {sites_str}")
            
            # Add concentration if available
            if concentration is not None:
                if concentration == 0:
                    output_lines.append(f"Concentration: 0.00e0 {self.polymer_data['concentration_units'] or ''}")
                else:
                    output_lines.append(f"Concentration: {concentration:.2e} {self.polymer_data['concentration_units'] or ''}")
            
            # Add free energy if available
            if free_energy is not None:
                output_lines.append(f"Free energy: {free_energy}")
            
            output_lines.append("")  # Empty line between polymers
        
        return '\n'.join(output_lines)