"""Module for reading and writing .tbnpolys files.

The .tbnpolys format provides a user-friendly text representation of polymers.
Format specification:
- Comments are designated by "#"
- Empty lines separate different polymers
- Each polymer is represented by its monomers, one per line
- Monomers can have a multiplicity prefix "n | " where n is the count
- Monomers are specified by name or binding site representation
"""

import re
from pathlib import Path
from typing import Any, List, Optional, Tuple

from .model import TBN, Monomer


class TbnpolysParser:
    """Parser for .tbnpolys files."""

    def __init__(self, tbn: Optional[TBN] = None):
        """Initialize parser with optional TBN context.

        Args:
            tbn: Optional TBN model for resolving monomer names and binding sites
        """
        self.tbn = tbn

    def parse_file(self, file_path: Path) -> List[List[Tuple[int, Any]]]:
        """Parse a .tbnpolys file.

        Args:
            file_path: Path to the .tbnpolys file

        Returns:
            List of polymers, where each polymer is a list of (multiplicity, monomer) tuples.
            The monomer can be either a string (name or binding sites) or a Monomer object if TBN is available.
        """
        with open(file_path) as f:
            content = f.read()
        return self.parse_content(content)

    def parse_content(self, content: str) -> List[List[Tuple[int, Any]]]:
        """Parse .tbnpolys content.

        Args:
            content: Content of a .tbnpolys file

        Returns:
            List of polymers, where each polymer is a list of (multiplicity, monomer) tuples.
        """
        polymers = []
        current_polymer = []

        lines = content.split("\n")
        i = 0
        while i < len(lines):
            original_line = lines[i]
            line = original_line.strip()

            # Check if this is a comment-only line (shouldn't count as empty)
            is_comment_only = line.startswith("#")

            # Remove comments from non-comment-only lines
            if not is_comment_only and "#" in line:
                line = line[: line.index("#")].strip()

            # Empty line signals end of current polymer (but not comment-only lines)
            if not line and not is_comment_only:
                if current_polymer:
                    polymers.append(current_polymer)
                    current_polymer = []
            elif not is_comment_only and line:
                # Parse monomer line
                multiplicity, monomer_spec = self._parse_monomer_line(line)
                if self.tbn:
                    monomer = self._resolve_monomer(monomer_spec)
                    current_polymer.append((multiplicity, monomer))
                else:
                    current_polymer.append((multiplicity, monomer_spec))

            i += 1

        # Add the last polymer if exists
        if current_polymer:
            polymers.append(current_polymer)

        return polymers

    def _parse_monomer_line(self, line: str) -> Tuple[int, str]:
        """Parse a single monomer line.

        Args:
            line: Line containing monomer specification

        Returns:
            Tuple of (multiplicity, monomer_spec)
        """
        # Check for multiplicity prefix "n | "
        match = re.match(r"^(\d+)\s*\|\s*(.+)$", line)
        if match:
            multiplicity = int(match.group(1))
            monomer_spec = match.group(2).strip()
        else:
            multiplicity = 1
            monomer_spec = line

        return multiplicity, monomer_spec

    def _resolve_monomer(self, monomer_spec: str) -> Monomer:
        """Resolve a monomer specification to a Monomer object.

        Args:
            monomer_spec: Either a monomer name, binding site representation,
                         or "name: binding_sites" format

        Returns:
            Monomer object

        Raises:
            ValueError: If monomer cannot be resolved or if name and binding sites don't match
        """
        if not self.tbn:
            raise ValueError("TBN context required to resolve monomers")

        # Check if monomer_spec has name: binding_sites syntax
        if ":" in monomer_spec:
            parts = monomer_spec.split(":", 1)
            if len(parts) == 2:
                name = parts[0].strip()
                binding_sites_str = parts[1].strip()

                # Try to find monomer by name
                for monomer in self.tbn.monomers:
                    if monomer.name == name:
                        # Verify that binding sites match
                        provided_sites = sorted(binding_sites_str.split())
                        monomer_sites = []
                        for site in monomer.binding_sites:
                            monomer_sites.append(site.name + ("*" if site.is_star else ""))
                        monomer_sites = sorted(monomer_sites)

                        if provided_sites != monomer_sites:
                            raise ValueError(
                                f"Monomer '{name}' exists but binding sites don't match. "
                                f"Expected: {' '.join(monomer_sites)}, "
                                f"Got: {' '.join(provided_sites)}"
                            )
                        return monomer

                # Name not found, raise error
                raise ValueError(f"Monomer with name '{name}' not found in TBN file")

        # First check if it's a monomer name (without colon syntax)
        for monomer in self.tbn.monomers:
            if monomer.name == monomer_spec:
                return monomer

        # Try to parse as binding sites
        binding_sites = monomer_spec.split()
        for monomer in self.tbn.monomers:
            # Check if binding sites match (order doesn't matter according to spec)
            monomer_sites = []
            for site in monomer.binding_sites:
                monomer_sites.append(site.name + ("*" if site.is_star else ""))

            if sorted(binding_sites) == sorted(monomer_sites):
                return monomer

        raise ValueError(f"Cannot resolve monomer: {monomer_spec}")


class TbnpolysWriter:
    """Writer for .tbnpolys files."""

    def __init__(self, tbn: TBN):
        """Initialize writer with TBN context.

        Args:
            tbn: TBN model for monomer information
        """
        self.tbn = tbn

    def write_polymers(
        self,
        polymers: List[List[int]],
        file_path: Path,
        concentrations: Optional[List[float]] = None,
        units: Optional[str] = None,
        header_comment: Optional[str] = None,
    ):
        """Write polymers to a .tbnpolys file.

        Args:
            polymers: List of polymer vectors (monomer counts)
            file_path: Output file path
            concentrations: Optional list of polymer concentrations
            units: Optional concentration units
            header_comment: Optional header comment
        """
        content = self.format_polymers(polymers, concentrations, units, header_comment)
        with open(file_path, "w") as f:
            f.write(content)

    def format_polymers(
        self,
        polymers: List[List[int]],
        concentrations: Optional[List[float]] = None,
        units: Optional[str] = None,
        header_comment: Optional[str] = None,
    ) -> str:
        """Format polymers as .tbnpolys content.

        Args:
            polymers: List of polymer vectors (monomer counts)
            concentrations: Optional list of polymer concentrations
            units: Optional concentration units
            header_comment: Optional header comment

        Returns:
            Formatted .tbnpolys content
        """
        lines = []

        # Add header comment if provided
        if header_comment:
            for line in header_comment.split("\n"):
                lines.append(f"# {line}")
            lines.append("")

        # Format each polymer
        for i, polymer in enumerate(polymers):
            polymer_lines = self._format_single_polymer(polymer)

            # Add concentration as comment if provided
            if concentrations and i < len(concentrations):
                conc_str = self._format_concentration(concentrations[i], units)
                polymer_lines.append(f"# Concentration: {conc_str}")

            lines.extend(polymer_lines)
            lines.append("")  # Empty line between polymers

        # Remove trailing empty line
        if lines and lines[-1] == "":
            lines.pop()

        return "\n".join(lines)

    def _format_single_polymer(self, polymer: List[int]) -> List[str]:
        """Format a single polymer.

        Args:
            polymer: Polymer vector (monomer counts)

        Returns:
            List of lines representing the polymer
        """
        lines = []

        for monomer_idx, count in enumerate(polymer):
            if count > 0:
                monomer = self.tbn.monomers[monomer_idx]
                monomer_spec = self._get_monomer_spec(monomer)

                if count == 1:
                    lines.append(monomer_spec)
                else:
                    lines.append(f"{count} | {monomer_spec}")

        return lines

    def _get_monomer_spec(self, monomer: Monomer) -> str:
        """Get the specification string for a monomer.

        Args:
            monomer: Monomer object

        Returns:
            Monomer name if available, otherwise binding site representation
        """
        if monomer.name:
            return monomer.name
        else:
            # Format binding sites
            sites = []
            for site in monomer.binding_sites:
                sites.append(site.name + ("*" if site.is_star else ""))
            return " ".join(sites)

    def _format_concentration(self, concentration: float, units: Optional[str]) -> str:
        """Format a concentration value nicely.

        Args:
            concentration: Concentration value
            units: Optional units string

        Returns:
            Formatted concentration string
        """
        # Format concentration nicely (avoid scientific notation for reasonable values)
        if concentration == 0:
            conc_str = "0"
        elif concentration >= 1000:
            conc_str = f"{concentration:.1e}"
        elif concentration >= 100:
            conc_str = f"{concentration:.1f}"
        elif concentration >= 10:
            conc_str = f"{concentration:.2f}"
        elif concentration >= 1:
            conc_str = f"{concentration:.3f}"
        elif concentration >= 0.01:
            conc_str = f"{concentration:.4f}"
        else:
            conc_str = f"{concentration:.2e}"

        if units:
            conc_str += f" {units}"

        return conc_str
