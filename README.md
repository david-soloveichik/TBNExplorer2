# TBN Explorer 2

A Python library and command-line tool for analyzing Thermodynamics of Binding Networks (TBN) models. This tool computes polymer bases (Hilbert bases) for TBN systems, which represent the fundamental "unsplittable" polymers in the network.

## Installation

### Prerequisites

1. **Python 3.8+** is required
2. **Normaliz** - A tool for discrete convex geometry computations
   - Download from: https://github.com/Normaliz/Normaliz
   - Default expected path: `/Users/dsolov/Documents/ResearchTools/Normaliz/normaliz`
   - You can specify a custom path using the `--normaliz-path` option

### Install from source

```bash
# Clone or download the repository
cd TBNExplorer2-CLI

# Install in development mode
pip install -e .

# Or install normally
pip install .
```

### Install dependencies only

```bash
pip install -r requirements.txt
```

## Usage

### Basic Command

```bash
tbnexplorer2 input.tbn
```

This will:
1. Parse the TBN file
2. Check the star-limiting restriction
3. Compute the polymer basis using Normaliz
4. Save results to `input-polymer-basis.txt`

### Command-Line Options

```bash
tbnexplorer2 input.tbn [options]

Options:
  -h, --help            Show help message
  -o, --output FILE     Specify output file (default: [input]-polymer-basis.txt)
  --normaliz-path PATH  Path to Normaliz executable
  -v, --verbose         Enable verbose output
  --check-only          Only check star-limiting restriction, don't compute basis
```

### Examples

```bash
# Basic usage
tbnexplorer2 examples/balanced.tbn

# Specify custom output file
tbnexplorer2 examples/balanced.tbn --output results.txt

# Verbose mode to see detailed progress
tbnexplorer2 examples/balanced.tbn --verbose

# Only check if TBN satisfies star-limiting restriction
tbnexplorer2 examples/balanced.tbn --check-only

# Use custom Normaliz installation
tbnexplorer2 examples/balanced.tbn --normaliz-path /usr/local/bin/normaliz
```

## TBN File Format

TBN files (`.tbn`) describe monomers and their binding sites. Here's the format:

### Basic Syntax

```
# Comments start with #
monomer_name: site1 site2 site3   # Named monomer
site1 site2 site3                  # Unnamed monomer
site1 site2 site3, 100.5          # Monomer with concentration
```

### Binding Sites

- **Unstar sites**: Regular binding sites (e.g., `a`, `b1`, `xyz`)
- **Star sites**: Complementary binding sites ending with `*` (e.g., `a*`, `b1*`, `xyz*`)
- Star and unstar sites with the same base name can bind together

### Rules

1. **Consistency**: Either ALL monomers have concentrations or NONE do
2. **Star-limiting**: For each binding site type, total unstar ≥ total star
3. **No special characters**: Binding sites cannot contain `,`, `*`, `|`, or `:`

### Example TBN Files

**Simple balanced TBN:**
```
# balanced.tbn - Each binding site has equal star/unstar
A: a b c
B: a* b* c*
C: d d e
D: d* d* e*
```

**TBN with concentrations:**
```
# with_concentrations.tbn
MA: a b, 100
MB: a* b*, 50
MC: c d, 75
MD: c* d*, 60
```

## Output Format

The polymer basis is saved in a human-readable format:

```
# Polymer basis - N polymers
#
# Polymer 1
2 | monomer_name
1 | a b c

# Polymer 2
1 | monomer_name
...
```

Each line shows:
- Count (`n |`) - multiplicity of the monomer in the polymer
- Monomer representation - either its name or binding sites

## Python API

You can also use TBN Explorer as a Python library:

```python
from tbnexplorer2 import TBNParser, TBN, PolymerBasisComputer

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
polymers = computer.compute_polymer_basis()

# Save results
computer.save_polymer_basis(polymers, "output.txt")
print(f"Found {len(polymers)} polymers in basis")
```

## Testing

Run the test suite:

```bash
# Install pytest if needed
pip install pytest

# Run all tests
pytest tests/

# Run with verbose output
pytest tests/ -v

# Run specific test file
pytest tests/test_parser.py
```

## Troubleshooting

### "Normaliz not found" Error

If you get an error about Normaliz not being found:

1. Install Normaliz from https://github.com/Normaliz/Normaliz
2. Either:
   - Use `--normaliz-path` to specify the location
   - Edit `NORMALIZ_PATH` in `tbnexplorer2/normaliz.py`

### Star-limiting Restriction Errors

If your TBN fails the star-limiting check:
- Ensure for each binding site type, you have at least as many unstar as star sites
- Use `--verbose` to see which binding sites are problematic
- Adjust concentrations if using them

### Installation Issues

If `pip install -e .` fails:
- Make sure you have Python 3.8 or newer
- Try `python -m pip install --upgrade pip` first
- Install numpy manually: `pip install numpy`

## License

This project is part of the TBN research toolkit for analyzing thermodynamics of binding networks.