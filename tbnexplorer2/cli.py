#!/usr/bin/env python3
"""
TBN Explorer 2 - Command Line Interface

A tool for analyzing Thermodynamics of Binding Networks (TBN) models.
"""

import argparse
import sys
import os
from pathlib import Path
from .parser import TBNParser
from .model import TBN
from .polymer_basis import PolymerBasisComputer
from .normaliz import NormalizRunner, NORMALIZ_PATH


def main():
    """Main entry point for the CLI."""
    parser = argparse.ArgumentParser(
        description="TBN Explorer 2 - Analyze Thermodynamics of Binding Networks",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Examples:
  # Compute polymer basis for a TBN file
  tbnexplorer2 example.tbn
  
  # Specify custom output file
  tbnexplorer2 example.tbn --output my-polymer-basis.txt
  
  # Use custom Normaliz path
  tbnexplorer2 example.tbn --normaliz-path /path/to/normaliz
  
  # Verbose output
  tbnexplorer2 example.tbn --verbose
        """
    )
    
    parser.add_argument(
        'input_file',
        help='Input TBN file'
    )
    
    parser.add_argument(
        '--output', '-o',
        help='Output file for polymer basis (default: [input]-polymer-basis.txt)'
    )
    
    parser.add_argument(
        '--normaliz-path',
        default=NORMALIZ_PATH,
        help=f'Path to Normaliz executable (default: {NORMALIZ_PATH})'
    )
    
    parser.add_argument(
        '--verbose', '-v',
        action='store_true',
        help='Enable verbose output'
    )
    
    parser.add_argument(
        '--check-only',
        action='store_true',
        help='Only check star-limiting restriction, do not compute polymer basis'
    )
    
    args = parser.parse_args()
    
    # Validate input file
    if not os.path.exists(args.input_file):
        print(f"Error: Input file '{args.input_file}' not found", file=sys.stderr)
        sys.exit(1)
    
    # Determine output file name
    if args.output:
        output_file = args.output
    else:
        base_name = Path(args.input_file).stem
        output_file = f"{base_name}-polymer-basis.txt"
    
    try:
        # Parse TBN file
        if args.verbose:
            print(f"Parsing TBN file: {args.input_file}")
        
        monomers, binding_site_index = TBNParser.parse_file(args.input_file)
        
        if args.verbose:
            print(f"Found {len(monomers)} monomers")
            print(f"Found {len(binding_site_index)} unique binding sites")
        
        # Create TBN model
        tbn = TBN(monomers, binding_site_index)
        
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
        
        # Check Normaliz availability
        normaliz_runner = NormalizRunner(args.normaliz_path)
        
        if not normaliz_runner.check_normaliz_available():
            print(f"Error: Normaliz not found at '{args.normaliz_path}'", file=sys.stderr)
            print("Please install Normaliz or specify the correct path with --normaliz-path", file=sys.stderr)
            sys.exit(1)
        
        # Compute polymer basis
        if args.verbose:
            print("Computing polymer basis...")
            print(f"Matrix A shape: {tbn.matrix_A.shape}")
        
        computer = PolymerBasisComputer(tbn, normaliz_runner)
        polymers = computer.compute_polymer_basis()
        
        if args.verbose:
            print(f"Found {len(polymers)} polymers in the basis")
        
        # Save polymer basis
        computer.save_polymer_basis(polymers, output_file)
        
        # Print summary
        print(f"Polymer basis computation complete")
        print(f"Number of polymers in basis: {len(polymers)}")
        print(f"Results saved to: {output_file}")
        
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