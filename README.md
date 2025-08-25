# TBN Explorer 2

A Python library and command-line tools for analyzing Thermodynamics of Binding Networks (TBN) models.

## Overview

The TBN (Thermodynamics of Binding Networks) model is a framework for studying abstract chemical systems at equilibrium. In this model:

- **Binding sites** (domains) are the fundamental units that can bind to complementary partners
- **Monomers** are collections (multisets) of binding sites  
- **Polymers** (complexes) are collections of monomers that bind together through complementary binding sites

Enthalpy is assumed to be infinitely preferred over entropy and the system must be "star-binding site limiting". Please see [1,2] for details of the TBN model.

Given a set of monomers (described by their binding sites) and their initial concentrations, TBN Explorer computes:
1. The **polymer basis** - the set of fundamental "unsplittable" polymers [3]
2. The **free energies** of each polymer equal to the number of bonds
3. The **equilibrium concentrations** of all polymers in the system [4]

### Extensions

TBN Explorer also includes advanced extensions for:
- **Canonical Reactions Analysis**: Identifying irreducible reactions between on-target and off-target polymers [5]
- **IBOT Algorithm**: Computing concentration exponents for maintaining detailed balance in polymer networks [5]

### References
[1] David Doty, Trent A. Rogers, David Soloveichik, Chris Thachuk, Damien Woods, "Thermodynamic binding networks," DNA23, 2017.
[2] Keenan Breik, Chris Thachuk, Marijn Heule, David Soloveichik, "Computing properties of stable configurations of thermodynamic binding networks," Theoretical Computer Science 785, 17-29, 2019.
[3] David Haley, David Doty, "Computing properties of thermodynamic binding networks: An integer programming approach," DNA27, 2021.
[4] https://coffeesolver.dev/
[5] Hamidreza Akef, Minki Hhan, David Soloveichik, "Computing and Bounding Equilibrium Concentrations in Athermic Chemical Systems," DNA31, 2025.


## Installation

### Prerequisites

1. **Python 3.8+** is required
2. **Python dependencies** - Automatically installed with pip:
   - numpy, scipy - Numerical computations
   - simpleeval - Safe expression evaluation for parametrized .tbn files
3. **Normaliz** - Tool for computing Hilbert bases and discrete convex geometry
   - Download from: https://github.com/Normaliz/Normaliz
   - Configure path via `NORMALIZ_PATH` environment variable or `.env` file
   - You can also specify a custom path using the `--normaliz-path` option
4. **COFFEE** - Tool for computing chemical equilibrium concentrations
   - Configure path via `COFFEE_CLI_PATH` environment variable or `.env` file
   - Required for equilibrium concentration calculations
5. **4ti2** (optional) - Alternative tool for computing Hilbert bases
   - Download from: https://4ti2.github.io/
   - Configure path via `FOURTI2_PATH` environment variable or `.env` file
   - Use with `--use-4ti2` flag as an alternative to Normaliz

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
  --use-4ti2                        Use 4ti2 instead of Normaliz for Hilbert basis computation
  --4ti2-path PATH                  Path to 4ti2 installation directory
  --user-friendly-polymer-basis     Generate human-readable polymer basis file
  --no-concentrations               Skip concentration calculations
  --no-free-energies                Skip free energy calculations (also skips concentrations)
  --parametrized var1=val1 ...      Provide values for template variables in .tbn file
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
  --min-concentration CONC   Minimum concentration threshold
  --constraints-file FILE    Advanced filtering with constraint specifications
```

**Constraint File Format:**
```
CONTAINS monomer1 monomer2    # Polymers containing these monomers
EXACTLY monomer1 monomer2     # Polymers with exactly these monomers
```

### tbnexplorer2-ibot - Iterative Balancing of Off-Target Polymers

Analyzes canonical reactions and assigns concentration exponents to off-target polymers using the IBOT algorithm.

```bash
tbnexplorer2-ibot input.tbn on_target.tbnpolys [options]

Options:
  --normaliz-path PATH              Path to Normaliz executable
  --use-4ti2                        Use 4ti2 instead of Normaliz
  --4ti2-path PATH                  Path to 4ti2 installation directory
  --generate-tbn CONC UNITS         Generate .tbn file with computed concentrations
                                    (e.g., --generate-tbn 0.01 nM)
  -v, --verbose                     Enable verbose output
```

**Purpose:**
- Identifies irreducible canonical reactions between on-target and off-target polymers
- Assigns concentration exponents to off-target polymers maintaining detailed balance
- Generates concentration specifications for all polymers

**Inputs:**
- `input.tbn`: TBN file without concentrations (no \UNITS)
- `on_target.tbnpolys`: File specifying which polymers are considered "on-target"

**Outputs:**
- `input-ibot-output.tbnpolys`: All polymers with concentration exponents (μ values)
- `input-ibot-cX.tbn`: Generated TBN file with concentrations (if --generate-tbn used)

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

### Parametrized TBN Files

TBN files support template variables and arithmetic expressions for flexible concentration specifications:

```
\UNITS: nM
# Simple variable substitution
monomer1: a b, {{conc1}}

# Arithmetic expressions
scaled: c d, {{base_conc * scale_factor}}
average: e f, {{(conc1 + conc2) / 2}}
complex: g h, {{conc1 * 0.7 + conc2 * 0.3}}

# Mix templates with literal values
fixed: i j, 100
```

**Running with parameters:**
```bash
tbnexplorer2 input.tbn --parametrized conc1=50 conc2=100 base_conc=75 scale_factor=2
```

**Supported operations in templates:**
- Basic arithmetic: `+`, `-`, `*`, `/`
- Exponentiation: `**` (e.g., `{{2 ** 3}}` evaluates to 8)
- Parentheses for grouping: `{{(a + b) * c}}`
- Floating point: `{{x * 1.5 + y * 0.5}}`

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
initiator
3 | propagator
```
Note: when the multiplicity prefix "{n} |" is missing, the multiplicity is assumed to be 1.

## Extensions

TBN Explorer includes advanced analysis tools for studying reaction networks and polymer dynamics.

### Canonical Reactions Analysis

The canonical reactions module identifies and analyzes irreducible reactions that generate off-target polymers from on-target ones.

**Key Concepts:**
- **On-target polymers**: Desired polymer configurations specified in a .tbnpolys file
- **Off-target polymers**: All other polymers in the polymer basis
- **Canonical reactions**: Reactions that create off-target polymers from purely on-target reactants
- **Irreducible canonical reactions**: Fundamental reactions that cannot be decomposed into simpler ones

### IBOT Algorithm (Iterative Balancing of Off-Target Polymers)

The IBOT algorithm assigns concentration exponents to polymers to maintain detailed balance in the system.

**Algorithm Overview:**
1. All on-target polymers are assigned concentration exponent μ = 1
2. Iteratively assigns concentration exponents to off-target polymers based on:
   - **Novelty**: Number of unassigned off-target polymers in a reaction
   - **Imbalance**: Difference in concentration exponents between reactants and products
   - **Imbalance-novelty ratio**: Used to determine assignment priority
3. Polymers unreachable from on-target reactions are excluded

**Example Workflow:**
```bash
# 1. Compute polymer basis from TBN without concentrations
tbnexplorer2 system.tbn --user-friendly-polymer-basis

# 2. Create on_target.tbnpolys file specifying desired polymers
# (manually or programmatically select from polymer basis)

# 3. Run IBOT analysis
tbnexplorer2-ibot system.tbn on_target.tbnpolys --verbose

# 4. Generate TBN with computed concentrations
tbnexplorer2-ibot system.tbn on_target.tbnpolys --generate-tbn 0.01 nM

# 5. Verify equilibrium with generated concentrations
tbnexplorer2 system-ibot-c0.01.tbn

# 6. Filter and view results
tbnexplorer2-filter system-ibot-c0.01.tbn -n 50
```

### .tbnpolys File Format

Polymer specification files for defining on-target polymers and viewing results.

**Format:**
```
# Comment describing polymer
2 | monomer_name     # 2 copies of named monomer
1 | a b c            # 1 copy of monomer with sites a, b, c
                     # Empty line separates polymers
# Next polymer
1 | monomer2
3 | d e f

# Concentration exponents (in IBOT output)
# μ: 1.5
```

**Key Rules:**
- Empty lines (not comments) separate different polymers
- Multiplicity prefix `n |` indicates n copies of the monomer
- Monomers can be specified by name or binding site list
- μ values appear in IBOT output files

## Python API

Use TBN Explorer as a Python library:

```python
from tbnexplorer2 import TBNParser, TBN, PolymerBasisComputer, FreeEnergyCalculator

# Parse TBN file
monomers, binding_sites = TBNParser.parse_file("input.tbn")

# Parse TBN file with parametrized variables
variables = {"conc1": 50, "conc2": 100, "scale": 1.5}
monomers, binding_sites, units, used_vars = TBNParser.parse_file("parametrized.tbn", variables)

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

### Extensions API

```python
from extensions import CanonicalReactionsComputer, IBOTAlgorithm
from tbnexplorer2 import TBNPolysParser

# Load on-target polymers
on_target_polymers = TBNPolysParser.parse_file("on_target.tbnpolys", tbn)

# Compute irreducible canonical reactions
reactions_computer = CanonicalReactionsComputer(tbn, polymer_basis)
irreducible_reactions = reactions_computer.compute_irreducible_canonical_reactions(
    on_target_polymers
)

# Run IBOT algorithm
ibot = IBOTAlgorithm(tbn, polymer_basis, irreducible_reactions)
concentration_exponents = ibot.compute_concentration_exponents(on_target_polymers)

# Generate monomer concentrations
monomer_concentrations = ibot.generate_monomer_concentrations(
    concentration_exponents, 
    base_concentration=0.01,  # in Molar
    units="nM"
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

# Test IBOT pipeline end-to-end
python extensions/test_ibot_pipeline.py

# Check code quality
ruff check .
ruff format .
```

### Testing the IBOT Pipeline

The `test_ibot_pipeline.py` script validates the complete IBOT workflow:
1. Generates concentration exponents using IBOT algorithm
2. Creates a TBN file with computed concentrations
3. Runs equilibrium calculation with tbnexplorer2
4. Verifies that equilibrium concentrations match expected values

This ensures the mathematical consistency of the concentration exponent assignments.

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