#!/usr/bin/env python3
"""
Custom argcomplete completers for TBNExplorer2 CLI tools.

This module provides specialized completion functions for arguments
that need custom logic beyond simple choices.
"""

import re
from pathlib import Path
from typing import Any, List

try:
    from argcomplete.completers import FilesCompleter
except ImportError:
    # Fallback if argcomplete is not installed
    FilesCompleter = None


# Built-in file completers for common file types
if FilesCompleter:
    TBNFilesCompleter = FilesCompleter(allowednames=(".tbn",), directories=True)
    TBNPolysFilesCompleter = FilesCompleter(allowednames=(".tbnpolys",), directories=True)
    TBNPolymatFilesCompleter = FilesCompleter(allowednames=(".tbnpolymat",), directories=True)
    TextFilesCompleter = FilesCompleter(allowednames=(".txt",), directories=True)
else:
    # Fallback to None if argcomplete is not installed
    TBNFilesCompleter = None
    TBNPolysFilesCompleter = None
    TBNPolymatFilesCompleter = None
    TextFilesCompleter = None


def concentration_units_completer(prefix: str = "", **kwargs) -> List[str]:
    """Complete concentration units."""
    units = ["nM", "pM", "uM", "mM", "M"]
    return [u for u in units if u.startswith(prefix)]


def parametrized_completer(prefix: str = "", parsed_args: Any = None, **kwargs) -> List[str]:
    """Complete parametrized variable assignments."""
    # Provide example format hints
    if not prefix or "=" not in prefix:
        return ["VAR=VALUE", "conc1=90.3", "conc2=50", "base_conc=100"]
    return []


def monomer_names_completer(prefix: str = "", parsed_args: Any = None, **kwargs) -> List[str]:
    """Extract and complete monomer names from the specified .tbn file."""
    if not parsed_args or not hasattr(parsed_args, "tbn_file"):
        return []

    tbn_file = getattr(parsed_args, "tbn_file", None)
    if not tbn_file:
        return []

    tbn_path = Path(tbn_file)
    if not tbn_path.exists() or not tbn_path.is_file():
        return []

    monomer_names = set()
    try:
        with open(tbn_path) as f:
            for line in f:
                # Skip comments
                line = line.split("#")[0].strip()
                if not line:
                    continue

                # Look for monomer names (format: "name:" or ">name")
                # Match "name:" at the beginning
                match = re.match(r"^\s*([^:>\s]+):", line)
                if match:
                    monomer_names.add(match.group(1))

                # Match ">name" anywhere in the line
                match = re.search(r">([^,\s]+)", line)
                if match:
                    monomer_names.add(match.group(1))
    except Exception:
        # If we can't read the file, just return empty list
        return []

    # Filter by prefix and return sorted list
    return sorted([name for name in monomer_names if name.startswith(prefix)])


def normaliz_path_completer(prefix: str = "", **kwargs) -> List[str]:
    """Complete paths to Normaliz executable."""

    # If user is typing a path, use file completion
    if "/" in prefix:
        # Let the default file completer handle this
        return []

    # Otherwise suggest common executable names
    common_names = ["normaliz", "Normaliz", "normaliz.exe"]
    return [name for name in common_names if name.startswith(prefix)]


def fourtitwo_path_completer(prefix: str = "", **kwargs) -> List[str]:
    """Complete paths to 4ti2 installation directory."""
    # This should complete directories
    # Let the default directory completer handle this
    return []


def coffee_path_completer(prefix: str = "", **kwargs) -> List[str]:
    """Complete paths to COFFEE executable."""

    # If user is typing a path, use file completion
    if "/" in prefix:
        # Let the default file completer handle this
        return []

    # Otherwise suggest common executable names
    common_names = ["coffee-cli", "coffee", "COFFEE"]
    return [name for name in common_names if name.startswith(prefix)]


def nupack_path_completer(prefix: str = "", **kwargs) -> List[str]:
    """Complete paths to NUPACK concentrations executable."""

    # If user is typing a path, use file completion
    if "/" in prefix:
        # Let the default file completer handle this
        return []

    # Otherwise suggest common executable names
    common_names = ["concentrations", "nupack-concentrations"]
    return [name for name in common_names if name.startswith(prefix)]
