"""
Filter module for TBN Explorer 2.

Filters polymers from .tbnpolymat files based on monomer name criteria.
"""

from pathlib import Path
from typing import List, Optional, Tuple

import numpy as np

from .parser import TBNParser
from .polymat_io import PolymatData, PolymatReader


def format_concentration_nicely(value: float, units: str) -> str:
    """
    Format concentration value nicely without scientific notation where possible.

    Args:
        value: Concentration value
        units: Concentration units (e.g., 'nM', 'uM', 'mM', 'M', 'pM')

    Returns:
        Nicely formatted string like "99.9 nM" instead of "9.99e+01 nM"
    """
    if value == 0:
        return f"0.00 {units}"

    # Use abs value for formatting, preserve sign separately
    abs_value = abs(value)
    sign = "-" if value < 0 else ""

    # For very small values (< 0.001), use scientific notation
    if abs_value < 0.001:
        return f"{sign}{abs_value:.2e} {units}"
    # For values between 0.001 and 0.01, show 3 significant figures after decimal
    elif abs_value < 0.01:
        return f"{sign}{abs_value:.4f} {units}"
    # For values between 0.01 and 0.1, show 3 significant figures after decimal
    elif abs_value < 0.1:
        return f"{sign}{abs_value:.3f} {units}"
    # For values between 0.1 and 1, show 2 decimal places
    elif abs_value < 1 or abs_value < 10:
        return f"{sign}{abs_value:.2f} {units}"
    # For values between 10 and 100, show 1 decimal place
    elif abs_value < 100 or abs_value < 1000:
        return f"{sign}{abs_value:.1f} {units}"
    # For values between 1000 and 10000, show no decimal places
    elif abs_value < 10000:
        return f"{sign}{abs_value:.0f} {units}"
    # For very large values, use scientific notation
    else:
        return f"{sign}{abs_value:.2e} {units}"


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
        self.monomers, self.binding_site_index, self.units = TBNParser.parse_file(str(self.tbn_file))

        # Validate that UNITS keyword exists (required for tbnexplorer2-filter)
        if self.units is None:
            raise ValueError(
                f"tbnexplorer2-filter requires a .tbn file with UNITS keyword and concentrations. "
                f"File '{self.tbn_file}' does not have UNITS specified."
            )

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

    def _load_polymat_file(self) -> PolymatData:
        """
        Load and parse the .tbnpolymat file.

        Returns:
            PolymatData object containing polymer information
        """
        reader = PolymatReader(str(self.polymat_file))
        return reader.read()

    def filter_by_monomers(
        self, monomer_names: List[str], percent_limit: Optional[float] = None, max_count: Optional[int] = None
    ) -> List[Tuple[int, np.ndarray, Optional[float], Optional[float]]]:
        """
        Filter polymers containing all specified monomers.

        Args:
            monomer_names: List of monomer names to filter by (can have duplicates for multiplicity).
                          If empty, returns all polymers (subject to other limits).
                          Note: Only named monomers can be filtered; unnamed monomers cannot be specified.
            percent_limit: Optional percentage limit (0-100) for filtering by concentration
            max_count: Maximum number of polymers to return

        Returns:
            List of tuples (polymer_index, monomer_counts, free_energy, concentration)
            sorted by decreasing concentration (if available) or by polymer index
        """
        # Handle empty monomer names list - return all polymers
        if not monomer_names:
            required_counts = {}
            monomer_name_to_indices = {}
        else:
            # Count required multiplicity for each monomer name
            required_counts = {}
            for name in monomer_names:
                required_counts[name] = required_counts.get(name, 0) + 1

            # Create mapping from monomer name to indices (ONLY for named monomers)
            monomer_name_to_indices = {}
            for i, monomer in enumerate(self.monomers):
                # Only use monomer name, not binding sites
                if monomer.name:
                    if monomer.name not in monomer_name_to_indices:
                        monomer_name_to_indices[monomer.name] = []
                    monomer_name_to_indices[monomer.name].append(i)

            # Check if all required monomer names exist
            for name in required_counts:
                if name not in monomer_name_to_indices:
                    # Return empty list if monomer name doesn't exist
                    return []

        # Calculate total concentration if available
        total_concentration = None
        if self.polymer_data.has_concentrations and self.polymer_data.concentrations is not None:
            total_concentration = np.sum(self.polymer_data.concentrations)

        # Filter polymers
        matching_polymers = []
        for i, polymer_counts in enumerate(self.polymer_data.polymers):
            # Check if polymer contains all required monomers with correct multiplicity
            matches = True
            if required_counts:  # Only check if we have filtering criteria
                for monomer_name, required_count in required_counts.items():
                    # Sum counts for all monomers with this name
                    actual_count = sum(polymer_counts[idx] for idx in monomer_name_to_indices[monomer_name])
                    if actual_count < required_count:
                        matches = False
                        break

            if matches:
                _, free_energy, concentration = self.polymer_data.get_polymer_data(i)

                # Apply percent limit if specified
                if percent_limit is not None and concentration is not None and total_concentration is not None:
                    concentration_percent = (concentration / total_concentration) * 100
                    if concentration_percent < percent_limit:
                        continue

                matching_polymers.append((i, polymer_counts, free_energy, concentration))

        # Sort by concentration (descending) if available
        if self.polymer_data.has_concentrations:
            matching_polymers.sort(key=lambda x: x[3] if x[3] is not None else 0, reverse=True)

        # Apply max_count limit
        if max_count is not None and max_count > 0:
            matching_polymers = matching_polymers[:max_count]

        return matching_polymers

    def _parse_constraints_file(self, filename: str) -> List[Tuple[str, List[str]]]:
        """
        Parse a constraints file.

        Args:
            filename: Path to the constraints file

        Returns:
            List of (constraint_type, monomer_names) tuples

        Raises:
            ValueError: If file format is invalid
        """
        import shlex

        constraints = []

        with open(filename) as f:
            for line_num, line in enumerate(f, 1):
                # Strip whitespace and skip empty lines or comments
                line = line.strip()
                if not line or line.startswith("#"):
                    continue

                # Parse constraint type and monomer names
                # Use shlex to handle quoted strings properly
                parts = shlex.split(line)
                if len(parts) < 1:
                    continue

                constraint_type = parts[0].upper()
                if constraint_type not in ["CONTAINS", "EXACTLY"]:
                    raise ValueError(
                        f"Invalid constraint type '{parts[0]}' on line {line_num}. Must be CONTAINS or EXACTLY."
                    )

                monomer_names = parts[1:]
                constraints.append((constraint_type, monomer_names))

        return constraints

    def _check_polymer_matches_constraint(
        self, polymer_counts: np.ndarray, constraint_type: str, required_monomer_names: List[str]
    ) -> bool:
        """
        Check if a polymer matches a single constraint.

        Args:
            polymer_counts: Array of monomer counts in the polymer
            constraint_type: 'CONTAINS' or 'EXACTLY'
            required_monomer_names: List of monomer names required by the constraint

        Returns:
            True if polymer matches the constraint
        """
        # Count required multiplicity for each monomer name
        required_counts = {}
        for name in required_monomer_names:
            required_counts[name] = required_counts.get(name, 0) + 1

        # Create mapping from monomer name to indices (ONLY for named monomers)
        monomer_name_to_indices = {}
        for i, monomer in enumerate(self.monomers):
            # Only use monomer name, not binding sites
            if monomer.name:
                if monomer.name not in monomer_name_to_indices:
                    monomer_name_to_indices[monomer.name] = []
                monomer_name_to_indices[monomer.name].append(i)

        # Check if all required monomer names exist
        for name in required_counts:
            if name not in monomer_name_to_indices:
                return False

        if constraint_type == "CONTAINS":
            # Check if polymer contains all required monomers with correct multiplicity
            for monomer_name, required_count in required_counts.items():
                actual_count = sum(polymer_counts[idx] for idx in monomer_name_to_indices[monomer_name])
                if actual_count < required_count:
                    return False
            return True

        elif constraint_type == "EXACTLY":
            # First check that all required monomers are present with exact multiplicity
            for monomer_name, required_count in required_counts.items():
                actual_count = sum(polymer_counts[idx] for idx in monomer_name_to_indices[monomer_name])
                if actual_count != required_count:
                    return False

            # Then check that no other monomers are present
            accounted_indices = set()
            for monomer_name in required_counts:
                accounted_indices.update(monomer_name_to_indices[monomer_name])

            return all(not (i not in accounted_indices and count > 0) for i, count in enumerate(polymer_counts))

        return False

    def filter_by_constraints_file(
        self, constraints_file: str, percent_limit: Optional[float] = None, max_count: Optional[int] = None
    ) -> List[Tuple[int, np.ndarray, Optional[float], Optional[float]]]:
        """
        Filter polymers based on constraints from a file.

        Args:
            constraints_file: Path to file containing constraints
            percent_limit: Optional percentage limit (0-100) for filtering by concentration
            max_count: Maximum number of polymers to return

        Returns:
            List of tuples (polymer_index, monomer_counts, free_energy, concentration)
            sorted by decreasing concentration
        """
        # Parse constraints from file
        constraints = self._parse_constraints_file(constraints_file)

        # If no constraints, return all polymers (subject to other limits)
        if not constraints:
            return self.filter_by_monomers([], percent_limit=percent_limit, max_count=max_count)

        # Calculate total concentration if available
        total_concentration = None
        if self.polymer_data.has_concentrations and self.polymer_data.concentrations is not None:
            total_concentration = np.sum(self.polymer_data.concentrations)

        # Filter polymers using OR logic for multiple constraints
        matching_polymers = []
        for i, polymer_counts in enumerate(self.polymer_data.polymers):
            # Check if polymer matches ANY constraint (OR logic)
            matches_any = False
            for constraint_type, monomer_names in constraints:
                if self._check_polymer_matches_constraint(polymer_counts, constraint_type, monomer_names):
                    matches_any = True
                    break

            if matches_any:
                _, free_energy, concentration = self.polymer_data.get_polymer_data(i)

                # Apply percent limit if specified
                if percent_limit is not None and concentration is not None and total_concentration is not None:
                    concentration_percent = (concentration / total_concentration) * 100
                    if concentration_percent < percent_limit:
                        continue

                matching_polymers.append((i, polymer_counts, free_energy, concentration))

        # Sort by concentration (descending) if available
        if self.polymer_data.has_concentrations:
            matching_polymers.sort(key=lambda x: x[3] if x[3] is not None else 0, reverse=True)

        # Apply max_count limit
        if max_count is not None and max_count > 0:
            matching_polymers = matching_polymers[:max_count]

        return matching_polymers

    def format_output_with_constraints(
        self,
        filtered_polymers: List[Tuple[int, np.ndarray, Optional[float], Optional[float]]],
        constraints_description: str,
        percent_limit: Optional[float] = None,
        max_count: Optional[int] = None,
    ) -> str:
        """
        Format filtered polymers for user-friendly output when using constraints file.

        Args:
            filtered_polymers: List of filtered polymer tuples
            constraints_description: Description of constraints used
            percent_limit: Percent limit used for filtering (if any)
            max_count: Maximum count limit used for filtering (if any)

        Returns:
            Formatted string output
        """
        output_lines = []

        # Header
        output_lines.append(f"# Filtered polymers using {constraints_description}")
        if percent_limit is not None:
            output_lines.append(f"# Percent limit: {percent_limit}%")
        if max_count is not None:
            output_lines.append(f"# Maximum count limit: {max_count}")
        output_lines.append(f"# Number of matching polymers: {len(filtered_polymers)}")

        # Calculate total concentration stats if available
        if self.polymer_data.has_concentrations and self.polymer_data.concentrations is not None:
            total_concentration = np.sum(self.polymer_data.concentrations)
            matching_concentration = sum(p[3] for p in filtered_polymers if p[3] is not None)
            percentage = (matching_concentration / total_concentration * 100) if total_concentration > 0 else 0
            output_lines.append(f"# Total concentration fraction: {percentage:.2f}%")
            if self.units:
                output_lines.append(f"# Concentration units: {self.units}")

        output_lines.append("#")

        # Format each polymer
        for idx, (_, monomer_counts, _, concentration) in enumerate(filtered_polymers, 1):
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

            # Add concentration if available (using nice formatting)
            if concentration is not None:
                units = self.units or ""
                formatted_conc = format_concentration_nicely(concentration, units)
                output_lines.append(f"Concentration: {formatted_conc}")

            # Note: Free energy is not included in output to save space

            output_lines.append("")  # Empty line between polymers

        return "\n".join(output_lines)

    def format_output(
        self,
        filtered_polymers: List[Tuple[int, np.ndarray, Optional[float], Optional[float]]],
        monomer_names: List[str],
        percent_limit: Optional[float] = None,
        max_count: Optional[int] = None,
    ) -> str:
        """
        Format filtered polymers for user-friendly output.

        Args:
            filtered_polymers: List of filtered polymer tuples
            monomer_names: Original list of monomer names used for filtering
            percent_limit: Percent limit used for filtering (if any)
            max_count: Maximum count limit used for filtering (if any)

        Returns:
            Formatted string output
        """
        output_lines = []

        # Header
        if monomer_names:
            output_lines.append(f"# Filtered polymers containing: {' '.join(monomer_names)}")
        else:
            output_lines.append("# All polymers")
        if percent_limit is not None:
            output_lines.append(f"# Percent limit: {percent_limit}%")
        if max_count is not None:
            output_lines.append(f"# Maximum count limit: {max_count}")
        output_lines.append(f"# Number of matching polymers: {len(filtered_polymers)}")

        # Calculate total concentration stats if available
        if self.polymer_data.has_concentrations and self.polymer_data.concentrations is not None:
            total_concentration = np.sum(self.polymer_data.concentrations)
            matching_concentration = sum(p[3] for p in filtered_polymers if p[3] is not None)
            percentage = (matching_concentration / total_concentration * 100) if total_concentration > 0 else 0
            output_lines.append(f"# Total concentration fraction: {percentage:.2f}%")
            if self.units:
                output_lines.append(f"# Concentration units: {self.units}")

        output_lines.append("#")

        # Format each polymer
        for idx, (_, monomer_counts, _, concentration) in enumerate(filtered_polymers, 1):
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

            # Add concentration if available (using nice formatting)
            if concentration is not None:
                units = self.units or ""
                formatted_conc = format_concentration_nicely(concentration, units)
                output_lines.append(f"Concentration: {formatted_conc}")

            # Note: Free energy is not included in output to save space

            output_lines.append("")  # Empty line between polymers

        return "\n".join(output_lines)
