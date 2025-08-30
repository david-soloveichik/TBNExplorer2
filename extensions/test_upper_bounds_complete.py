#!/usr/bin/env python3
"""Test that upper bounds computation with all off-target polymers gives identical results to regular IBOT."""

import re
import sys
from pathlib import Path


def parse_tbnpolys_with_mu(filepath):
    """Parse .tbnpolys file and extract polymer descriptions and their μ values."""
    polymers = {}
    current_polymer = []
    current_mu = None

    with open(filepath) as f:
        for line in f:
            line = line.strip()

            # Skip comments except μ lines
            if line.startswith("#"):
                if line.startswith("# μ:"):
                    # Extract μ value
                    mu_match = re.match(r"# μ: ([-\d.]+)", line)
                    if mu_match:
                        current_mu = float(mu_match.group(1))
                continue

            # Empty line indicates end of polymer
            if not line:
                if current_polymer and current_mu is not None:
                    # Sort monomers to create canonical representation
                    polymer_key = tuple(sorted(current_polymer))
                    polymers[polymer_key] = current_mu
                    current_polymer = []
                    current_mu = None
            else:
                # Add monomer to current polymer
                current_polymer.append(line)

    # Handle last polymer if file doesn't end with empty line
    if current_polymer and current_mu is not None:
        polymer_key = tuple(sorted(current_polymer))
        polymers[polymer_key] = current_mu

    return polymers


def main():
    """Compare results from regular IBOT vs upper bounds with all off-target polymers."""

    # File paths
    regular_file = "and_gate_noA-ibot.tbnpolys"
    upper_bounds_file = "and_gate_noA-ibot-upper-bounds.tbnpolys"

    # Check files exist
    if not Path(regular_file).exists():
        print(f"Error: {regular_file} not found")
        sys.exit(1)

    if not Path(upper_bounds_file).exists():
        print(f"Error: {upper_bounds_file} not found")
        sys.exit(1)

    # Parse both files
    print(f"Parsing {regular_file}...")
    regular_polymers = parse_tbnpolys_with_mu(regular_file)

    print(f"Parsing {upper_bounds_file}...")
    upper_bounds_polymers = parse_tbnpolys_with_mu(upper_bounds_file)

    print(f"\nRegular IBOT found {len(regular_polymers)} polymers with μ values")
    print(f"Upper bounds IBOT found {len(upper_bounds_polymers)} polymers with μ values")

    # Compare polymer sets
    regular_keys = set(regular_polymers.keys())
    upper_bounds_keys = set(upper_bounds_polymers.keys())

    only_in_regular = regular_keys - upper_bounds_keys
    only_in_upper_bounds = upper_bounds_keys - regular_keys
    common_polymers = regular_keys & upper_bounds_keys

    if only_in_regular:
        print(f"\nPolymers only in regular IBOT: {len(only_in_regular)}")
        for polymer in sorted(only_in_regular):
            print(f"  {polymer}: μ = {regular_polymers[polymer]}")

    if only_in_upper_bounds:
        print(f"\nPolymers only in upper bounds IBOT: {len(only_in_upper_bounds)}")
        for polymer in sorted(only_in_upper_bounds):
            print(f"  {polymer}: μ = {upper_bounds_polymers[polymer]}")

    # Compare μ values for common polymers
    print(f"\nComparing μ values for {len(common_polymers)} common polymers:")

    all_identical = True
    tolerance = 1e-10

    for polymer in sorted(common_polymers):
        mu_regular = regular_polymers[polymer]
        mu_upper = upper_bounds_polymers[polymer]
        diff = abs(mu_regular - mu_upper)

        if diff > tolerance:
            all_identical = False
            print(f"  DIFFERENCE: {polymer}")
            print(f"    Regular: μ = {mu_regular}")
            print(f"    Upper bounds: μ = {mu_upper}")
            print(f"    Difference: {diff}")

    if all_identical and len(regular_polymers) == len(upper_bounds_polymers):
        print("\n✅ SUCCESS: All polymers have identical μ values!")
        print("This confirms the upper bounds implementation is correct.")
    else:
        print("\n❌ FAILURE: Results differ between regular and upper bounds IBOT")
        if len(regular_polymers) != len(upper_bounds_polymers):
            print(f"  Different number of polymers: {len(regular_polymers)} vs {len(upper_bounds_polymers)}")

    # Print summary of μ values
    print("\n=== Summary of μ values ===")

    # Collect unique μ values from regular run
    unique_mu_regular = sorted(set(regular_polymers.values()))
    print(f"\nRegular IBOT unique μ values: {unique_mu_regular}")

    # Collect unique μ values from upper bounds run
    unique_mu_upper = sorted(set(upper_bounds_polymers.values()))
    print(f"Upper bounds IBOT unique μ values: {unique_mu_upper}")

    return 0 if all_identical else 1


if __name__ == "__main__":
    sys.exit(main())
