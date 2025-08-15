# TBN Explorer 2

A Python library and command-line tools for analyzing Thermodynamics of Binding Networks (TBN) models.

## Overview

The TBN (Thermodynamics of Binding Networks) model is a framework for studying abstract chemical systems at equilibrium. In this model:

- **Binding sites** (domains) are the fundamental units that can bind to complementary partners
- **Monomers** are collections (multisets) of binding sites  
- **Polymers** (complexes) are collections of monomers that bind together through complementary binding sites

The system must be star limiting and enthalpy is assumed to be infinite.

Given a set of monomers (described by their binding sites) and their initial concentrations, TBN Explorer computes:
1. The **polymer basis** - the set of fundamental "unsplittable" polymers
2. The **free energies** of each polymer equal to the number of bonds
3. The **equilibrium concentrations** of all polymers in the system



## Installation

### Prerequisites

1. **Python 3.8+** is required
2. **Normaliz** - Tool for computing Hilbert bases and discrete convex geometry
   - Download from: https://github.com/Normaliz/Normaliz
   - Configure path via `NORMALIZ_PATH` environment variable or `.env` file
   - You can also specify a custom path using the `--normaliz-path` option
3. **COFFEE** - Tool for computing chemical equilibrium concentrations
   - Configure path via `COFFEE_CLI_PATH` environment variable or `.env` file
   - Required for equilibrium concentration calculations

### Install from source

```bash
# Clone or download the repository
cd TBNExplorer2-CLI

# Copy and configure environment file
cp .env.example .env
# Edit .env to set NORMALIZ_PATH and COFFEE_CLI_PATH

# Install in development mode
pip install -e .

# Or install normally
pip install .
```

## Command-Line Tools

### tbnexplorer2 - Main Analysis Tool

Computes polymer bases, free energies, and equilibrium concentrations.

```bash
tbnexplorer2 input.tbn [options]

Options:
  -h, --help                        Show help message
  --normaliz-path PATH              Path to Normaliz executable
  --user-friendly-polymer-basis     Generate human-readable polymer basis file
  --no-concentrations               Skip concentration calculations
  --no-free-energies                Skip free energy calculations (also skips concentrations)
  -v, --verbose                     Enable verbose output
```

**Outputs:**
- `input.tbnpolymat` - Matrix file with polymer compositions, free energies, and concentrations
- `input-polymer-basis.txt` - Human-readable polymer basis (with `--user-friendly-polymer-basis`)

### tbnexplorer2-filter - Polymer Filtering Tool

Filters and displays polymers containing specific monomers.

```bash
tbnexplorer2-filter input.tbn [monomer1] [monomer2] ... [options]

Options:
  -n, --num N                Maximum number of polymers to output (default: 100)
  -p, --percent-limit P      Only show polymers above P% of total concentration
  --constraints-file FILE    Advanced filtering with constraint specifications
```

**Constraint File Format:**
```
CONTAINS monomer1 monomer2    # Polymers containing these monomers
EXACTLY monomer1 monomer2     # Polymers with exactly these monomers
```

## TBN File Format (.tbn)

TBN files describe monomers as collections of binding sites with optional concentrations.

### Basic Syntax

```
# Comments start with #
\UNITS: nM                         # Concentration units (nM, pM, uM, mM, M)

monomer_name: site1 site2 site3, 100   # Named monomer with concentration
site1 site2* >another_name, 50         # Alternative naming syntax
site1 site2 site3                      # Unnamed monomer (no concentration without UNITS)
```

### Binding Sites

- **Unstar sites**: Regular binding sites (e.g., `a`, `b1`, `xyz`)
- **Star sites**: Complementary sites ending with `*` (e.g., `a*`, `b1*`)
- Complementary pairs (`a` and `a*`) can bind together
- Prohibited characters in names: `,`, `>`, `*`, `|`, `:`, `\`

### Key Rules

1. **Units Consistency**: With `\UNITS` keyword, ALL monomers must have concentrations. Without it, NO monomers can have concentrations.
2. **Star-limiting Restriction**: For each binding site type, total unstar count ≥ total star count across all monomers
3. **Duplicate Monomers**: 
   - Without UNITS: Treated independently
   - With UNITS: Concentrations are summed (must have same name or one named)

### Example TBN Files

**Simple system without concentrations:**
```
# Simple binding network
A: a b c
B: a* b* c*
dimer: d d e
complement: d* d* e*
```

**System with concentrations:**
```
\UNITS: nM
initiator: a b, 100
propagator: a* b c, 50
terminator: c* d, 75
cap: d*, 25.5
```

## Output Formats

### .tbnpolymat File

Matrix format with polymers sorted by concentration:
```
# Header with metadata
\MATRIX-HASH: [hash]
\UNITS: nM
# Columns: [monomer multiplicities] [free energy] [concentration]
1 0 2 0 -3 45.2
0 1 1 1 -2 12.8
...
```

### User-Friendly Polymer Basis

Human-readable format showing polymer composition:
```
# Polymer 1
2 | monomer_name
1 | a b c

# Polymer 2
1 | initiator
3 | propagator
```

## Python API

Use TBN Explorer as a Python library:

```python
from tbnexplorer2 import TBNParser, TBN, PolymerBasisComputer, FreeEnergyCalculator

# Parse TBN file
monomers, binding_sites = TBNParser.parse_file("input.tbn")

# Create TBN model
tbn = TBN(monomers, binding_sites)

# Check star-limiting restriction
is_valid, error_msg = tbn.check_star_limiting()
if not is_valid:
    print(f"Error: {error_msg}")
    exit(1)

# Compute polymer basis
computer = PolymerBasisComputer(tbn)
polymer_basis = computer.compute_polymer_basis()

# Calculate free energies
calc = FreeEnergyCalculator(tbn)
free_energies = calc.compute_free_energies(polymer_basis)

# Compute equilibrium concentrations (requires COFFEE)
if tbn.has_concentrations:
    concentrations = computer.compute_equilibrium_concentrations(
        polymer_basis, free_energies, tbn.monomer_concentrations
    )
```

## Performance Considerations

### Caching

The polymer basis computation (most expensive operation) is cached based on the monomer matrix hash. When re-running with different concentrations but same monomers, the cached basis is reused automatically.

### Large Systems

TBN Explorer is designed to handle:
- Polymer bases with hundreds of thousands of polymers
- Efficient matrix operations using NumPy
- Streaming I/O for large output files

## Testing

Run the test suite:

```bash
# Run all tests
pytest tests/

# Run with verbose output
pytest tests/ -v

# Run specific test
pytest tests/test_parser.py::TestTBNParser::test_units_parsing

# Check code quality
ruff check .
ruff format .
```

## Troubleshooting

### "Normaliz not found" Error

1. Install Normaliz from https://github.com/Normaliz/Normaliz
2. Either:
   - Set `NORMALIZ_PATH` in your `.env` file
   - Use `--normaliz-path` to specify location
   - Export `NORMALIZ_PATH` environment variable

### "COFFEE not found" Error

1. Build COFFEE from source (see COFFEE documentation)
2. Set `COFFEE_CLI_PATH` in your `.env` file or environment

### Star-limiting Restriction Errors

- Ensure for each binding site: unstar count ≥ star count
- Use `--verbose` to identify problematic sites
- Check concentration values if using UNITS

### Missing Concentrations

- With `\UNITS`: ALL monomers need concentrations
- Without `\UNITS`: NO monomers should have concentrations
- Use `--no-concentrations` to skip equilibrium calculations


## License

MIT License

Copyright (c) 2024 TBN Explorer Contributors

Permission is hereby granted, free of charge, to any person obtaining a copy
of this software and associated documentation files (the "Software"), to deal
in the Software without restriction, including without limitation the rights
to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
copies of the Software, and to permit persons to whom the Software is
furnished to do so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in all
copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
SOFTWARE.