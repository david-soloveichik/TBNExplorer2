#!/usr/bin/env python3
"""
TBN Explorer 2 Filter - Command Line Interface

Filters polymers from .tbnpolymat files based on monomer name criteria.
"""

import argparse
import sys
from pathlib import Path
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
  
  # Filter by unnamed monomers using their binding sites
  tbnexplorer2-filter example.tbn "a1* a2* b1* b2*"
  tbnexplorer2-filter example.tbn "a2 b1" "c1* c2*"
        """
    )
    
    parser.add_argument(
        'tbn_file',
        help='Input TBN file (corresponding .tbnpolymat file must exist)'
    )
    
    parser.add_argument(
        'monomer_names',
        nargs='+',
        help='Space-separated list of monomer names to filter by. For named monomers, use their name (e.g., B, C). For unnamed monomers, use their binding sites (e.g., "a1* a2* b1* b2*"). Duplicates increase required multiplicity.'
    )
    
    parser.add_argument(
        '--percent-limit',
        type=float,
        metavar='P',
        help='Only show polymers with concentration > P%% of total concentration'
    )
    
    args = parser.parse_args()
    
    # Validate input file
    if not Path(args.tbn_file).exists():
        print(f"Error: Input file '{args.tbn_file}' not found", file=sys.stderr)
        sys.exit(1)
    
    # Validate percent limit if provided
    if args.percent_limit is not None:
        if args.percent_limit < 0 or args.percent_limit > 100:
            print(f"Error: --percent-limit must be between 0 and 100", file=sys.stderr)
            sys.exit(1)
    
    try:
        # Create filter and load data
        polymer_filter = PolymerFilter(args.tbn_file)
        
        # Filter polymers
        filtered_polymers = polymer_filter.filter_by_monomers(
            args.monomer_names,
            percent_limit=args.percent_limit
        )
        
        # Format and output results
        output = polymer_filter.format_output(
            filtered_polymers,
            args.monomer_names,
            percent_limit=args.percent_limit
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