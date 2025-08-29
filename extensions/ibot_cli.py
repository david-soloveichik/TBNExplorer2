#!/usr/bin/env python3
"""Command-line interface for the IBOT algorithm.

This module provides the tbnexplorer2-ibot command for running the
Iterative Balancing of Off-Target (IBOT) algorithm.
"""

import argparse
import sys
from pathlib import Path

import numpy as np

from tbnexplorer2.fourtitwo import FourTiTwoRunner
from tbnexplorer2.model import TBN
from tbnexplorer2.normaliz import NormalizRunner
from tbnexplorer2.parser import TBNParser
from tbnexplorer2.polymer_basis import PolymerBasisComputer

from .canonical_reactions import CanonicalReactionsComputer
from .ibot import IBOTAlgorithm


def main():
    """Main entry point for tbnexplorer2-ibot CLI."""
    parser = argparse.ArgumentParser(description="Run IBOT algorithm for iterative balancing of off-target polymers")

    # Required arguments
    parser.add_argument("tbn_file", type=str, help="Input .tbn file (without concentrations)")
    parser.add_argument("on_target_file", type=str, help="Input .tbnpolys file specifying on-target polymers")

    # Optional arguments
    parser.add_argument(
        "--use-4ti2",
        action="store_true",
        help="Use 4ti2 instead of Normaliz for both Hilbert basis computations (polymer basis and irreducible canonical reactions)",
    )
    parser.add_argument(
        "--generate-tbn",
        nargs=2,
        metavar=("C", "UNITS"),
        help="Generate .tbn file with concentrations using base value C and specified units",
    )
    parser.add_argument(
        "--output-prefix", type=str, help="Prefix for output files (default: input filename without extension)"
    )
    parser.add_argument(
        "--output-canonical-reactions",
        action="store_true",
        help="Generate text file showing all irreducible canonical reactions ordered by IBOT iteration",
    )
    parser.add_argument(
        "--upper-bound-on-polymers",
        type=str,
        metavar="TBNPOLYS_FILE",
        help="Compute upper bounds only for specific off-target polymers listed in this .tbnpolys file. "
        "Could be much faster than full computation for large systems. "
        "Incompatible with --generate-tbn.",
    )

    args = parser.parse_args()

    # Validate input files
    tbn_path = Path(args.tbn_file)
    on_target_path = Path(args.on_target_file)

    if not tbn_path.exists():
        print(f"Error: TBN file not found: {tbn_path}", file=sys.stderr)
        sys.exit(1)

    if not on_target_path.exists():
        print(f"Error: On-target polymers file not found: {on_target_path}", file=sys.stderr)
        sys.exit(1)

    # Validate upper bounds options
    if args.upper_bound_on_polymers:
        upper_bound_path = Path(args.upper_bound_on_polymers)
        if not upper_bound_path.exists():
            print(f"Error: Upper bound polymers file not found: {upper_bound_path}", file=sys.stderr)
            sys.exit(1)

        if args.generate_tbn:
            print(
                "Error: --upper-bound-on-polymers cannot be used with --generate-tbn (we don't know all polymer concentrations)",
                file=sys.stderr,
            )
            sys.exit(1)

    # Determine output prefix
    output_prefix = args.output_prefix or tbn_path.stem

    try:
        # Step 1: Parse TBN file
        print(f"Parsing TBN file: {tbn_path}")
        monomers, binding_site_index, concentration_units, used_variables = TBNParser.parse_file(str(tbn_path))

        # Verify no concentrations are specified
        if concentration_units is not None:
            print("Error: TBN file must not contain concentrations (no \\UNITS)", file=sys.stderr)
            sys.exit(1)

        # Create TBN object
        tbn = TBN(monomers, binding_site_index, concentration_units)

        # Step 2: Compute polymer basis
        print("Computing polymer basis...")

        # Choose Hilbert basis solver
        if args.use_4ti2:
            runner = FourTiTwoRunner()
            print("Using 4ti2 for Hilbert basis computation")
        else:
            runner = NormalizRunner()
            print("Using Normaliz for Hilbert basis computation")

        basis_computer = PolymerBasisComputer(tbn, runner)
        polymers = basis_computer.compute_polymer_basis()
        polymer_vectors = [p.monomer_counts for p in polymers]
        print(f"Found {len(polymers)} polymers in the basis")

        # Step 3: Set up canonical reactions computation
        print("Setting up canonical reactions computation...")
        reactions_computer = CanonicalReactionsComputer(tbn, use_4ti2=args.use_4ti2)

        # Load on-target polymers
        on_target_indices = reactions_computer.load_on_target_polymers(on_target_path, polymer_vectors)
        print(f"Loaded {len(on_target_indices)} on-target polymers")

        # Set up matrices
        reactions_computer.setup_matrices(polymer_vectors, on_target_indices)

        # Step 4: Compute irreducible canonical reactions
        if args.upper_bound_on_polymers:
            # Load target polymers for upper bounds
            from tbnexplorer2.tbnpolys_io import TbnpolysParser

            print(f"Loading target polymers from {args.upper_bound_on_polymers}...")
            parser = TbnpolysParser(tbn)
            target_polymers_raw = parser.parse_file(args.upper_bound_on_polymers)

            # Convert to polymer indices
            target_polymer_indices = set()
            for polymer_raw in target_polymers_raw:
                # Create monomer count vector
                counts = np.zeros(len(tbn.monomers), dtype=int)
                for multiplicity, monomer in polymer_raw:
                    monomer_idx = tbn.monomers.index(monomer)
                    counts[monomer_idx] += multiplicity

                # Find index in polymer basis
                found = False
                for i, polymer in enumerate(polymer_vectors):
                    if np.array_equal(counts, polymer):
                        target_polymer_indices.add(i)
                        found = True
                        break

                if not found:
                    print(f"Warning: Target polymer {counts} not found in polymer basis", file=sys.stderr)

            if not target_polymer_indices:
                print("Error: No valid target polymers found in polymer basis", file=sys.stderr)
                sys.exit(1)

            print(f"Computing irreducible canonical reactions for {len(target_polymer_indices)} target polymers...")
            reactions = reactions_computer.compute_irreducible_canonical_reactions_for_targets(target_polymer_indices)
            print(f"Found {len(reactions)} irreducible canonical reactions that produce target polymers")
        else:
            print("Computing irreducible canonical reactions...")
            reactions = reactions_computer.compute_irreducible_canonical_reactions()
            print(f"Found {len(reactions)} irreducible canonical reactions")

        # Step 5: Check detailed balance of on-target polymers
        print("Checking detailed balance of on-target polymers...")
        violating_reaction = reactions_computer.check_on_target_detailed_balance(reactions)

        if violating_reaction:
            print("Error: On-target polymers not in detailed balance", file=sys.stderr)
            print(f"Violating reaction: {violating_reaction}", file=sys.stderr)
            sys.exit(1)

        print("On-target polymers are in detailed balance")

        # Step 6: Run IBOT algorithm
        print("\nRunning IBOT algorithm...")
        ibot = IBOTAlgorithm(tbn, polymer_vectors, on_target_indices, reactions)
        concentration_exponents = ibot.run()

        # Step 7: Generate output files
        suffix = "-upper-bounds" if args.upper_bound_on_polymers else ""
        output_tbnpolys = Path(f"{output_prefix}-ibot{suffix}.tbnpolys")
        print(f"\nGenerating output .tbnpolys file: {output_tbnpolys}")
        ibot.generate_tbnpolys_output(output_tbnpolys)

        # Generate reactions output file if requested
        if args.output_canonical_reactions:
            output_reactions = Path(f"{output_prefix}-ibot{suffix}-reactions.txt")
            print(f"\nGenerating canonical reactions output file: {output_reactions}")
            ibot.generate_reactions_output(output_reactions)

        # Generate .tbn file if requested
        if args.generate_tbn:
            c = float(args.generate_tbn[0])
            units = args.generate_tbn[1]

            # Validate units
            valid_units = ["nM", "pM", "uM", "mM", "M"]
            if units not in valid_units:
                print(f"Error: Invalid units '{units}'. Must be one of: {', '.join(valid_units)}", file=sys.stderr)
                sys.exit(1)

            output_tbn = Path(f"{output_prefix}-ibot-c{c}.tbn")
            print(f"\nGenerating .tbn file with concentrations: {output_tbn}")
            ibot.generate_tbn_output(output_tbn, c, units)

        print("\nIBOT algorithm completed successfully!")

        # Print summary statistics
        n_on_target = len(on_target_indices)
        n_off_target = len(polymer_vectors) - n_on_target
        unique_mus = len(set(concentration_exponents.values()))

        print("\nSummary:")
        print(f"  Total polymers: {len(polymer_vectors)}")
        print(f"  On-target polymers: {n_on_target}")
        print(f"  Off-target polymers: {n_off_target}")
        print(f"  Unique concentration exponents: {unique_mus}")

    except Exception as e:
        print(f"Error: {e}", file=sys.stderr)
        import traceback

        traceback.print_exc()
        sys.exit(1)


if __name__ == "__main__":
    main()
