import re
from typing import List, Tuple, Optional, Dict
from .model import Monomer, BindingSite


class TBNParser:
    """Parser for TBN (Thermodynamics of Binding Networks) input files."""
    
    @staticmethod
    def parse_file(filepath: str) -> Tuple[List[Monomer], Dict[str, int]]:
        """
        Parse a TBN file and extract monomers with their concentrations.
        
        Args:
            filepath: Path to the TBN file
            
        Returns:
            Tuple of (list of Monomer objects, dict of binding site to index)
            
        Raises:
            ValueError: If file format is invalid or mixing concentration specifications
        """
        monomers = []
        binding_site_index = {}
        has_concentration = None
        line_number = 0
        
        with open(filepath, 'r') as f:
            for line in f:
                line_number += 1
                # Remove comments
                if '#' in line:
                    line = line[:line.index('#')]
                
                line = line.strip()
                
                # Skip empty lines
                if not line:
                    continue
                
                # Parse the line
                monomer_data = TBNParser._parse_monomer_line(line, line_number)
                if monomer_data:
                    name, binding_sites, concentration = monomer_data
                    
                    # Check consistency of concentration specification
                    if has_concentration is None:
                        has_concentration = (concentration is not None)
                    elif has_concentration != (concentration is not None):
                        raise ValueError(
                            f"Line {line_number}: Inconsistent concentration specification. "
                            "Either all monomers must have concentrations or none."
                        )
                    
                    # Update binding site index
                    for site in binding_sites:
                        base_site = site.name
                        if base_site not in binding_site_index:
                            binding_site_index[base_site] = len(binding_site_index)
                    
                    # Create monomer
                    monomer = Monomer(
                        name=name,
                        binding_sites=binding_sites,
                        concentration=concentration,
                        original_line=line
                    )
                    monomers.append(monomer)
        
        if not monomers:
            raise ValueError("No valid monomers found in file")
        
        return monomers, binding_site_index
    
    @staticmethod
    def _parse_monomer_line(line: str, line_number: int) -> Optional[Tuple[Optional[str], List[BindingSite], Optional[float]]]:
        """
        Parse a single monomer line.
        
        Format: [name:] site1 site2 ... [, concentration]
        
        Args:
            line: Line to parse
            line_number: Line number for error reporting
            
        Returns:
            Tuple of (name or None, list of BindingSite objects, concentration or None)
            None if line cannot be parsed
        """
        # Check for name (indicated by colon)
        name = None
        remaining = line
        
        if ':' in line:
            parts = line.split(':', 1)
            name = parts[0].strip()
            remaining = parts[1].strip()
            
            # Validate name doesn't contain special characters
            if any(c in name for c in ',*|:'):
                raise ValueError(f"Line {line_number}: Invalid monomer name '{name}' - cannot contain ,*|:")
        
        # Check for concentration (after comma)
        concentration = None
        if ',' in remaining:
            parts = remaining.rsplit(',', 1)
            remaining = parts[0].strip()
            try:
                concentration = float(parts[1].strip())
                if concentration < 0:
                    raise ValueError(f"Line {line_number}: Negative concentration not allowed")
            except ValueError as e:
                if "Negative concentration" in str(e):
                    raise
                raise ValueError(f"Line {line_number}: Invalid concentration value '{parts[1].strip()}'")
        
        # Parse binding sites
        site_strings = remaining.split()
        if not site_strings:
            return None
        
        binding_sites = []
        for site_str in site_strings:
            # Validate binding site string
            if any(c in site_str for c in ',|:'):
                raise ValueError(f"Line {line_number}: Invalid binding site '{site_str}' - cannot contain ,|:")
            
            # Check if it's a star binding site
            if site_str.endswith('*'):
                base_name = site_str[:-1]
                is_star = True
            else:
                base_name = site_str
                is_star = False
            
            if not base_name:
                raise ValueError(f"Line {line_number}: Invalid binding site '{site_str}'")
            
            binding_sites.append(BindingSite(base_name, is_star))
        
        return name, binding_sites, concentration