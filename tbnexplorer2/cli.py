#!/usr/bin/env python3
"""
TBN Explorer 2 - Command Line Interface

A tool for analyzing Thermodynamics of Binding Networks (TBN) models.
"""

import argparse
import os
import sys
from pathlib import Path

from .coffee import COFFEE_CLI_PATH, COFFEERunner
from .fourtitwo import FOURTI2_PATH, FourTiTwoRunner
from .model import TBN
from .normaliz import NORMALIZ_PATH, NormalizRunner
from .parser import TBNParser
from .polymer_basis import PolymerBasisComputer
from .units import get_unit_display_name


def main():
    """Main entry point for the CLI."""
    parser = argparse.ArgumentParser(
        description="TBN Explorer 2 - Analyze Thermodynamics of Binding Networks",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog=r"""
Examples:
  # Compute polymer basis for a TBN file (generates .tbnpolymat only)
  tbnexplorer2 example.tbn
  
  # Also save user-friendly polymer basis text file
  tbnexplorer2 example.tbn --user-friendly-polymer-basis
  
  # Specify custom output file for user-friendly basis
  tbnexplorer2 example.tbn --user-friendly-polymer-basis --output my-polymer-basis.txt
  
  # Use custom Normaliz path
  tbnexplorer2 example.tbn --normaliz-path /path/to/normaliz
  
  # Verbose output
  tbnexplorer2 example.tbn --verbose

TBN File Format:
  Concentration units are specified directly in the .tbn file using:
    \UNITS: nM    (or pM, uM, mM, M)
  
  If \UNITS is specified, ALL monomers must have concentrations.
  If \UNITS is not specified, NO monomers can have concentrations.
        """,
    )

    parser.add_argument("input_file", help="Input TBN file")

    parser.add_argument("--output", "-o", help="Output file for polymer basis (default: [input]-polymer-basis.txt)")

    parser.add_argument(
        "--normaliz-path", default=NORMALIZ_PATH, help=f"Path to Normaliz executable (default: {NORMALIZ_PATH})"
    )

    parser.add_argument(
        "--use-4ti2",
        action="store_true",
        help="Use 4ti2 instead of Normaliz for Hilbert basis computation",
    )

    parser.add_argument(
        "--4ti2-path", default=FOURTI2_PATH, help=f"Path to 4ti2 installation directory (default: {FOURTI2_PATH})"
    )

    parser.add_argument("--verbose", "-v", action="store_true", help="Enable verbose output")

    parser.add_argument(
        "--check-only", action="store_true", help="Only check star-limiting restriction, do not compute polymer basis"
    )

    parser.add_argument(
        "--no-concentrations",
        action="store_true",
        help="Do not compute equilibrium concentrations even if monomer concentrations are provided",
    )

    parser.add_argument(
        "--no-free-energies",
        action="store_true",
        help="Do not compute polymer free energies (also disables concentration computation)",
    )

    parser.add_argument(
        "--coffee-path", default=COFFEE_CLI_PATH, help=f"Path to COFFEE executable (default: {COFFEE_CLI_PATH})"
    )

    parser.add_argument(
        "--user-friendly-polymer-basis",
        action="store_true",
        help="Save user-friendly polymer basis to [input]-polymer-basis.txt file",
    )

    parser.add_argument(
        "--parametrized",
        nargs="*",
        metavar="VAR=VALUE",
        help="Variable assignments for parametrized .tbn files (e.g., a=20 b=100.4)",
    )

    args = parser.parse_args()

    # Validate input file
    if not os.path.exists(args.input_file):
        print(f"Error: Input file '{args.input_file}' not found", file=sys.stderr)
        sys.exit(1)

    # Determine output file names
    input_path = Path(args.input_file)
    base_name = input_path.stem

    output_file = args.output or str(input_path.parent / f"{base_name}-polymer-basis.txt")

    # Always generate .tbnpolymat file
    polymat_file = str(input_path.parent / f"{base_name}.tbnpolymat")

    # Parse parametrized arguments if provided
    variables = {}
    if args.parametrized:
        for assignment in args.parametrized:
            if "=" not in assignment:
                print(
                    f"Error: Invalid parameter assignment '{assignment}'. Expected format: VAR=VALUE", file=sys.stderr
                )
                sys.exit(1)
            var_name, var_value = assignment.split("=", 1)
            var_name = var_name.strip()
            try:
                variables[var_name] = float(var_value.strip())
            except ValueError:
                print(f"Error: Invalid numeric value '{var_value}' for parameter '{var_name}'", file=sys.stderr)
                sys.exit(1)

    try:
        # Parse TBN file
        if args.verbose:
            print(f"Parsing TBN file: {args.input_file}")
            if variables:
                print(f"Using parameters: {variables}")

        monomers, binding_site_index, concentration_units, used_variables = TBNParser.parse_file(
            args.input_file, variables=variables
        )

        if args.verbose:
            print(f"Found {len(monomers)} monomers")
            print(f"Found {len(binding_site_index)} unique binding sites")

        # Create TBN model
        tbn = TBN(monomers, binding_site_index, concentration_units=concentration_units)

        # Check star-limiting restriction
        if args.verbose:
            print("Checking star-limiting restriction...")

        is_valid, error_msg = tbn.check_star_limiting()

        if not is_valid:
            print(f"Error: {error_msg}", file=sys.stderr)
            sys.exit(1)

        if args.verbose:
            print("TBN satisfies star-limiting restriction")

        if args.check_only:
            print("Star-limiting check passed")
            sys.exit(0)

        # Choose Hilbert basis solver
        if args.use_4ti2:
            # Use 4ti2
            solver_runner = FourTiTwoRunner(getattr(args, "4ti2_path"))
            if not solver_runner.check_fourtitwo_available():
                print(f"Error: 4ti2 not found at '{getattr(args, '4ti2_path')}'", file=sys.stderr)
                print("Please install 4ti2 or specify the correct path with --4ti2-path", file=sys.stderr)
                sys.exit(1)
            if args.verbose:
                print("Using 4ti2 for Hilbert basis computation")
        else:
            # Use Normaliz (default)
            solver_runner = NormalizRunner(args.normaliz_path)
            if not solver_runner.check_normaliz_available():
                print(f"Error: Normaliz not found at '{args.normaliz_path}'", file=sys.stderr)
                print("Please install Normaliz or specify the correct path with --normaliz-path", file=sys.stderr)
                sys.exit(1)
            if args.verbose:
                print("Using Normaliz for Hilbert basis computation")

        # Try to load cached polymer basis first
        computer = PolymerBasisComputer(tbn, solver_runner)
        polymers = computer.load_cached_polymer_basis(polymat_file)

        if polymers is not None:
            used_cache = True
            print("Using cached polymer basis (matrix hashes match)")
            if args.verbose:
                print(f"Loaded {len(polymers)} polymers from cache")
        else:
            used_cache = False
            # Compute polymer basis
            print("Computing polymer basis...")
            if args.verbose:
                print(f"Matrix A shape: {tbn.matrix_A.shape}")

            polymers = computer.compute_polymer_basis()

            if args.verbose:
                print(f"Found {len(polymers)} polymers in the basis")

        # Save polymer basis in user-friendly format if requested
        if args.user_friendly_polymer_basis:
            computer.save_polymer_basis(polymers, output_file)

        # Determine computation options
        compute_free_energies = not args.no_free_energies
        compute_concentrations = not args.no_concentrations and compute_free_energies

        # Check if COFFEE is available if needed
        coffee_runner = None
        if compute_concentrations and tbn.concentrations is not None:
            coffee_runner = COFFEERunner(args.coffee_path)
            if not coffee_runner.check_coffee_available():
                print(f"Warning: COFFEE not found at '{args.coffee_path}', skipping concentration computation")
                compute_concentrations = False
                coffee_runner = None

        # Save .tbnpolymat file
        computer.save_tbnpolymat(
            polymers,
            polymat_file,
            compute_free_energies=compute_free_energies,
            compute_concentrations=compute_concentrations,
            coffee_runner=coffee_runner,
            verbose=args.verbose,
            parameters=used_variables,
        )

        # Print summary
        if used_cache:
            print("\nSummary:")
            print(f"Polymer basis: {len(polymers)} polymers (cached)")
        else:
            print("\nPolymer basis computation complete")
            print(f"Number of polymers in basis: {len(polymers)}")

        # Show concentration units information only if concentrations are provided
        if tbn.concentrations is not None:
            unit_display = get_unit_display_name(concentration_units)
            print(f"Concentration units: {unit_display}")

        print("Results saved to:")
        if args.user_friendly_polymer_basis:
            print(f"  - Polymer basis: {output_file}")
        print(f"  - Polymer matrix: {polymat_file}")

        if not compute_free_energies:
            print("Note: Free energies not computed (--no-free-energies flag)")
        elif tbn.concentrations is None:
            print("Note: No monomer concentrations provided, equilibrium concentrations not computed")
        elif not compute_concentrations:
            print("Note: Equilibrium concentrations not computed (--no-concentrations flag or COFFEE unavailable)")

    except ValueError as e:
        print(f"Error parsing TBN file: {e}", file=sys.stderr)
        sys.exit(1)
    except RuntimeError as e:
        print(f"Error during computation: {e}", file=sys.stderr)
        sys.exit(1)
    except Exception as e:
        print(f"Unexpected error: {e}", file=sys.stderr)
        if args.verbose:
            import traceback

            traceback.print_exc()
        sys.exit(1)


if __name__ == "__main__":
    main()
