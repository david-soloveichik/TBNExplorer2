import re
import numpy as np
from typing import List, Tuple, Optional, Dict
from .model import Monomer, BindingSite


class TBNParser:
    """Parser for TBN (Thermodynamics of Binding Networks) input files."""
    
    @staticmethod
    def parse_file(filepath: str) -> Tuple[List[Monomer], Dict[str, int], Optional[str]]:
        """
        Parse a TBN file and extract monomers with their concentrations.
        
        Args:
            filepath: Path to the TBN file
            
        Returns:
            Tuple of (list of Monomer objects, dict of binding site to index, concentration units or None)
            
        Raises:
            ValueError: If file format is invalid or UNITS/concentration specifications are inconsistent
        """
        # First pass: scan for UNITS keyword
        units = None
        with open(filepath, 'r') as f:
            for line_number, line in enumerate(f, 1):
                # Remove comments
                if '#' in line:
                    line = line[:line.index('#')]
                line = line.strip()
                
                # Skip empty lines
                if not line:
                    continue
                
                # Check for UNITS keyword
                if line.startswith('UNITS:'):
                    if units is not None:
                        raise ValueError(f"Line {line_number}: Multiple UNITS specifications found")
                    try:
                        units = line.split(':', 1)[1].strip()
                        if units not in ['nM', 'pM', 'uM', 'mM', 'M']:
                            raise ValueError(f"Line {line_number}: Invalid units '{units}'. Must be one of: nM, pM, uM, mM, M")
                    except IndexError:
                        raise ValueError(f"Line {line_number}: Invalid UNITS format. Expected 'UNITS: <unit>'")
                    continue
                
                # If we encounter a non-UNITS, non-comment line, break to start monomer parsing
                break
        
        # Second pass: parse monomers
        monomers = []
        binding_site_index = {}
        line_number = 0
        monomer_names = set()  # Track monomer names for uniqueness check
        
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
                
                # Skip UNITS lines during monomer parsing
                if line.startswith('UNITS:'):
                    continue
                
                # Parse the line
                monomer_data = TBNParser._parse_monomer_line(line, line_number, units)
                if monomer_data:
                    name, binding_sites, concentration = monomer_data
                    
                    # Check consistency of UNITS and concentration specifications
                    if units is not None:  # UNITS specified - all monomers must have concentrations
                        if concentration is None:
                            raise ValueError(
                                f"Line {line_number}: UNITS specified but monomer lacks concentration. "
                                "When UNITS is present, all monomers must have concentrations."
                            )
                    else:  # No UNITS - no monomers can have concentrations
                        if concentration is not None:
                            raise ValueError(
                                f"Line {line_number}: Monomer has concentration but no UNITS specified. "
                                "When concentrations are used, UNITS must be specified."
                            )
                    
                    # Check for conflicts between monomer names and binding sites
                    if name:
                        # Check if monomer name conflicts with existing binding sites
                        if name in binding_site_index:
                            raise ValueError(
                                f"Line {line_number}: Monomer name '{name}' conflicts with binding site name"
                            )
                        monomer_names.add(name)
                    
                    # Update binding site index and check for conflicts
                    for site in binding_sites:
                        base_site = site.name
                        # Check if binding site conflicts with existing monomer names
                        if base_site in monomer_names:
                            raise ValueError(
                                f"Line {line_number}: Binding site '{base_site}' conflicts with monomer name"
                            )
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
        
        # Handle monomer repetition when UNITS is present
        if units is not None:
            monomers = TBNParser._aggregate_identical_monomers(monomers, binding_site_index)
        
        return monomers, binding_site_index, units
    
    @staticmethod
    def _parse_monomer_line(line: str, line_number: int, units: Optional[str] = None) -> Optional[Tuple[Optional[str], List[BindingSite], Optional[float]]]:
        """
        Parse a single monomer line.
        
        Format: [name:] site1 site2 ... [, concentration]
        
        Args:
            line: Line to parse
            line_number: Line number for error reporting
            units: Units specification (allows negative concentrations if present)
            
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
            
            # Validate name doesn't contain special characters or spaces
            if any(c in name for c in ',*|:'):
                raise ValueError(f"Line {line_number}: Invalid monomer name '{name}' - cannot contain ,*|:")
            if ' ' in name:
                raise ValueError(f"Line {line_number}: Invalid monomer name '{name}' - cannot contain spaces")
        
        # Check for concentration (after comma)
        concentration = None
        if ',' in remaining:
            parts = remaining.rsplit(',', 1)
            remaining = parts[0].strip()
            try:
                concentration = float(parts[1].strip())
                # Only disallow negative concentrations when UNITS is not present
                if concentration < 0 and units is None:
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
    
    @staticmethod
    def _aggregate_identical_monomers(monomers: List[Monomer], binding_site_index: Dict[str, int]) -> List[Monomer]:
        """
        Aggregate identical monomers by summing their concentrations.
        
        Two monomers are considered identical if they have the same binding sites
        (order-independent comparison using vector representation).
        
        Args:
            monomers: List of monomers to aggregate
            binding_site_index: Dictionary mapping binding site names to indices
            
        Returns:
            List of aggregated monomers with summed concentrations
            
        Raises:
            ValueError: If any final aggregated concentration is negative
        """
        if not monomers:
            return monomers
            
        # Group monomers by their vector representation (which is order-independent)
        vector_to_monomers = {}
        
        for monomer in monomers:
            vector = monomer.to_vector(binding_site_index)
            vector_key = tuple(vector)  # Convert to tuple for use as dict key
            
            if vector_key not in vector_to_monomers:
                vector_to_monomers[vector_key] = []
            vector_to_monomers[vector_key].append(monomer)
        
        # Aggregate concentrations for identical monomers
        aggregated_monomers = []
        
        for vector_key, monomer_group in vector_to_monomers.items():
            if len(monomer_group) == 1:
                # No duplicates, keep as-is
                aggregated_monomers.append(monomer_group[0])
            else:
                # Sum concentrations of identical monomers
                total_concentration = sum(m.concentration for m in monomer_group)
                
                # Check for negative final concentration
                if total_concentration < 0:
                    raise ValueError(
                        f"Negative final concentration ({total_concentration}) after aggregating "
                        f"identical monomers: {monomer_group[0].get_binding_sites_str()}"
                    )
                
                # Create aggregated monomer using the first occurrence's properties
                # but with summed concentration
                first_monomer = monomer_group[0]
                aggregated_monomer = Monomer(
                    name=first_monomer.name,
                    binding_sites=first_monomer.binding_sites,
                    concentration=total_concentration,
                    original_line=f"# Aggregated from {len(monomer_group)} identical monomers"
                )
                aggregated_monomers.append(aggregated_monomer)
        
        return aggregated_monomers