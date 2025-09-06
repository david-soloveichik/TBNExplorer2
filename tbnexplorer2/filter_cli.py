#!/usr/bin/env python3
# PYTHON_ARGCOMPLETE_OK
"""
TBN Explorer 2 Filter - Command Line Interface

Filters polymers from .tbnpolymat files based on monomer name criteria.
"""

import argparse
import sys
from pathlib import Path

try:
    import argcomplete
except ImportError:
    argcomplete = None

from .completers import TBNFilesCompleter, TextFilesCompleter, monomer_names_completer
from .filter import PolymerFilter


def main():
    """Main entry point for the filter CLI."""
    parser = argparse.ArgumentParser(
        description="Filter polymers from .tbnpolymat files by monomer names",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Examples:
  # Filter polymers containing monomers B and C
  tbnexplorer2-filter example.tbn B C
  
  # Filter polymers with at least 2 copies of monomer B
  tbnexplorer2-filter example.tbn B B
  
  # Filter with concentration percent limit (show only polymers > 1.5% of total)
  tbnexplorer2-filter example.tbn B C --percent-limit 1.5
  
  # Filter using constraints file for advanced filtering
  tbnexplorer2-filter example.tbn --constraints-file constraints.txt
  
  # Show all polymers (no filtering)
  tbnexplorer2-filter example.tbn
        """,
    )

    tbn_arg = parser.add_argument("tbn_file", help="Input TBN file (corresponding .tbnpolymat file must exist)")
    if TBNFilesCompleter:
        tbn_arg.completer = TBNFilesCompleter

    monomer_arg = parser.add_argument(
        "monomer_names",
        nargs="*",
        help="Space-separated list of monomer names to filter by. Only named monomers can be filtered (e.g., B, C, output). Duplicates increase required multiplicity. If no monomers specified, returns all polymers (subject to other limits).",
    )
    if argcomplete:
        monomer_arg.completer = monomer_names_completer

    parser.add_argument(
        "--num", "-n", type=int, default=100, metavar="N", help="Maximum number of polymers to output (default: 100)"
    )

    parser.add_argument(
        "--percent-limit",
        "-p",
        type=float,
        metavar="P",
        help="Only show polymers with concentration > P%% of total concentration",
    )

    constraints_arg = parser.add_argument(
        "--constraints-file",
        type=str,
        metavar="FILE",
        help="File containing advanced filtering constraints (CONTAINS or EXACTLY)",
    )
    if TextFilesCompleter:
        constraints_arg.completer = TextFilesCompleter

    # Enable argcomplete if available
    if argcomplete:
        argcomplete.autocomplete(parser)

    args = parser.parse_args()

    # Validate input file
    if not Path(args.tbn_file).exists():
        print(f"Error: Input file '{args.tbn_file}' not found", file=sys.stderr)
        sys.exit(1)

    # Validate constraints file compatibility
    if args.constraints_file and args.monomer_names:
        print("Error: Cannot specify monomer names on command line when using --constraints-file", file=sys.stderr)
        sys.exit(1)

    # Validate constraints file exists if provided
    if args.constraints_file and not Path(args.constraints_file).exists():
        print(f"Error: Constraints file '{args.constraints_file}' not found", file=sys.stderr)
        sys.exit(1)

    # Validate percent limit if provided
    if args.percent_limit is not None and (args.percent_limit < 0 or args.percent_limit > 100):
        print("Error: --percent-limit must be between 0 and 100", file=sys.stderr)
        sys.exit(1)

    # Validate num parameter
    if args.num < 1:
        print("Error: --num must be at least 1", file=sys.stderr)
        sys.exit(1)

    try:
        # Create filter and load data
        polymer_filter = PolymerFilter(args.tbn_file)

        # Filter polymers
        if args.constraints_file:
            # Use constraints file filtering
            filtered_polymers = polymer_filter.filter_by_constraints_file(
                args.constraints_file, percent_limit=args.percent_limit, max_count=args.num
            )
            constraints_description = f"constraints from {Path(args.constraints_file).name}"
        else:
            # Use regular monomer name filtering
            filtered_polymers = polymer_filter.filter_by_monomers(
                args.monomer_names, percent_limit=args.percent_limit, max_count=args.num
            )
            constraints_description = None

        # Format and output results
        if args.constraints_file:
            output = polymer_filter.format_output_with_constraints(
                filtered_polymers, constraints_description, percent_limit=args.percent_limit, max_count=args.num
            )
        else:
            output = polymer_filter.format_output(
                filtered_polymers, args.monomer_names, percent_limit=args.percent_limit, max_count=args.num
            )

        print(output)

    except FileNotFoundError as e:
        print(f"Error: {e}", file=sys.stderr)
        sys.exit(1)
    except ValueError as e:
        print(f"Error: {e}", file=sys.stderr)
        sys.exit(1)
    except Exception as e:
        print(f"Unexpected error: {e}", file=sys.stderr)
        import traceback

        traceback.print_exc()
        sys.exit(1)


if __name__ == "__main__":
    main()
