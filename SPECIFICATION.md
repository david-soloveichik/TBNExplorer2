
# Background 
We are working with thermodynamics of abstract chemical systems (the model is called the TBN model). A TBN is defined in terms of binding sites (aka domains), monomers, and polymers (aka complexes). Conceptually, monomers are multisets of binding sites, and polymers are multisets of monomers. The big picture goal is given a set of monomers (described in terms of their binding sites) and their concentrations, we want to compute the equilibrium concentrations of the polymers. 

We will be building a Python library. The functionality should also be accessible via command line (CLI) tools including `tbnexplorer2` and other secondary tools.

# Input file encoding (.tbn) and TBN representation
A binding site is a symbol like "a" or "ahu34_8" or "5" (any string not containing white space nor any of the following: ">,*|:\\" which will have special meaning). There are two types of binding sites, unstar and star. Star binding sites always have a "*" as the last character. The idea is that corresponding star and unstar binding sites bind (they are complementary).
Examples: 
    Unstar binding sites: a, 2, ab, 12_a
    Corresponding (complementary) star binding sites: a*, 2*, 12_a*

Monomers are multisets of binding sites. For example: {a, b, a*, b}.

Internally, we represent a monomer as an integer vector. In this vector, we represent star binding sites as negative numbers. So for example, if we have four kinds of (unstarred) binding sites: (a, b, c, d) then the monomer {a, a*, d, b*, d} would be represented as: <0, -1, 0, 2>. Note that the "a" and "a*" effectively cancel in a single monomer because one adds the +1 and the other adds a -1. 

The set of monomers in our TBN corresponds to a matrix A of the monomer vectors. Specifically, each monomer will correspond to a *column* of the matrix A.

We need to be efficient in representing the matrix A as we will be doing some linear algebra on it (eg matrix multiplication and transpose).

The TBN is input from a .tbn file that is given as a command line argument. Here is an example of the input format. Empty lines, or lines just containing comments are ok.
```
# this line is a comment
\UNITS: nM
monomer1: a a* b2 b1, 100
b a c* >monomer2, 10.5
b2 b1* >C
b2* b2* b1 b1, 75.5  # this is a comment
a b2* b1, 50.7
```
Note that "#" followed by anything is a comment and should be ignored. 
The number following a comma is the concentration of that monomer. 
`<name>:` and `><name>` are two different ways to specify monomer names (see below).
`\UNITS:` is a keyword (see below). 

## Monomer names
Monomers can have names which is indicated by a name following by ":" prior to the monomer specification.
They can also be indicated by ">" followed by the monomer name after the monomer specification.
To avoid confusion, monomer names and binding sites should be distinct, and there should be an error otherwise. (Capitalization matters, so it's ok to have a binding site "c" and a monomer named "C", for example.) 
Monomer names should not have spaces or any of the prohibited symbols ">,*|:\\".

## Units and concentrations specification
The presence or absence of a `UNITS` keyword determines the type of .tbn file:
- **With UNITS**: All monomers MUST have concentrations specified. The UNITS line specifies the concentration units and has the format `\UNITS: <unit>` where `<unit>` can be: `nM` (nanoMolar), `pM` (picoMolar), `uM` (microMolar), `mM` (milliMolar), or `M` (Molar).
- **Without UNITS**: NO monomers can have concentrations specified.

The UNITS line can appear anywhere before the first monomer definition (comments and empty lines are allowed before it).

When UNITS is specified, monomers must have concentrations indicated by a _comma_ followed by the concentration after the monomer specification.
Important: An error should be returned if UNITS is present but some monomers lack concentrations, or if UNITS is absent but some monomers have concentrations.

## Monomer repetition:
The same monomer could occur multiple times in the .tbn file. 

- Without `\UNITS`: Treat the two identical monomers independently 

- With `\UNITS`: Have only one entry for the monomer in the matrix A. The monomer's concentration should be the _sum_ of the concentrations of its multiple entries. Important: this includes the possibility of a negative concentration of one of the entries. We should do error checking to make sure all the final concentrations of all monomers is non-negative. 
If there are duplicate monomers with different names, we should return an error. If one of the duplicate monomers has a name `<name>`, while the other duplicates are nameless, this is ok, and the resulting monomer with the summed concentration should have the name `<name>`.


# Polymers 
A polymer is a multiset of monomers. Note that we can have a polymer (or complex) consisting of just a single monomer (we would call this a "singleton polymer").

Internally, we will represent a polymer as vector of non-negative integers corresponding to the counts of the monomers in it. Given our representation of monomers in matrix A above, if x is a polymer, then b = A.x is a vector of binding site *excesses* in the polymer x. If b = 0, then every binding site is paired to its complement.


**Star-limiting restriction**: We always enforce the "star-limited restriction" described below, and give an error if it is not satisfied. Intuitively, a TBN is star-limited if for every binding site there is at least as much x as x* (totalled over all monomers). More formally:

- When monomer concentrations are given in the TBN file: Let c be the vector of the concentrations of the monomers. Then the TBN is *star-limited* if A.c ≥ 0 (component-wise).
- When monomer concentrations are not given, then c = Transpose[(1,1,...,1)] for the star-limited restriction check.


# Additional Utilities
- We want to have Python utility functions for running `Normaliz` to solve Hilbert basis problems and retrieving its output. It is a bit tricky to extract the Hilbert basis from the command line output of Normaliz.


# Computations over TBNs
Given an input TBN, we will be interested in computing several things. 

## 1. Polymer basis (aka Hilbert basis of polymers)
The polymer basis is conceptually the set of "unsplittable" polymers, the polymers that cannot be decomposed into two without losing some bonding. We compute it by solving a special linear problem. 

First, we need to make sure that for every binding site x, there is a column in A consisting of all zeros and a "-1" in the position corresponding to x. This is intuitively equivalent to having a singleton monomer {x*} in the TBN. If there are some binding sites for which such a singleton monomer does not exist in the TBN, then we add these columns to A as if we had more monomers. Importantly, we need to keep track of how many monomers we had originally, so that we can remove these fake added monomers later. Call this modified matrix A'. 

Now we want to find the Hilbert basis of the cone corresponding to solutions of: 
A'.x = 0
x ≥ 0 
We do this using `Normaliz`.  
Let H be the Hilbert basis output by Normaliz for this problem. To get the *polymer basis*, we remove the last entries of vectors in H corresponding to the "fake" singleton monomers (that we added to form A'), and remove any resulting duplicates in H. 

Important: The Hilbert basis H may be very large! There could be hundreds of thousands of vectors in H. Thus these operations need to be very efficient. 

We optionally save the polymer basis into a user-friendly text file. If the command line option `--user-friendly-polymer-basis` is given to out CLI tool `tbnexplorer2`, it should do the following: 
Assuming we are given `<example>.tbn` as input, it should create a file called `<example>-polymer-basis.tbnpolys`. This file should represent the polymers in the polymer basis in the user-friendly way as described below in the section on .tbnpolys files. Each polymer should start with a comment line like "# Polymer 3" to indicate that this is the 3rd polymer.

The CLI tool should return to the command line the number of polymers in the polymer basis.

### Alternative Hilbert basis calculator
Let's have an optional argument `--use-4ti2` to use the Hilbert basis solver of the `4ti2` tool instead of Normaliz. (We hope it might be faster in some situations.) 

## 2. Compute polymer free energies
Given the polymers in the polymer basis, we want to compute their "free energies". 
In general, we want to implement the following: Given a matrix of polymers, with a polymer per row (ie vector of the monomer counts of that polymer), we want to compute the "free energies" of each polymer. 
Since the binding is always maximized, we can set the free energy of a polymer x as dG(x) = 0.

However, sometimes we want to add an extra dG_assoc penalty term for each additional monomer in a polymer based on empirical parameters (as is done in Nupack). If the optional `--deltaG-assoc <dG_assoc> <dH_assoc>` parameter is given to `tbnexplorer2`, we increase (i.e., make less negative) the free energy of each polymer by `assoc_energy_penalty` as computed in `/PLANNING/assoc_energy.py`.

Variable mapping:
`total_monomers` = number of monomers in the polymer,
`temp_C` = temperature in Celcius (default 37C, can be changed by `--temp` described below)
`G_BIMOLECULAR` = `dG_assoc`,
`H_BIMOLECULAR` = `dH_assoc`.

## 3. Compute equilibrium polymer concentrations
The big picture is that we want to compute the equilibrium concentrations of all the polymers in the polymer basis. We use the command line tool COFFEE for this.

COFFEE takes two input files: CFE and CON. 
- The CFE file contains the following matrix: Each line corresponds to a polymer in the polymer basis (i.e., counts of monomers in that polymer), followed by the free energy of that polymer. 
- The CON file contains the concentration of each monomer, one per line. So there should be as many lines as columns in the CFE file minus one (the free energy). The monomer concentrations should be as given in the input .tbn file, _in units of Molar_ (see below).
Use the `-o <filename>` flag for `coffee-cli` to write its output in `<filename>`. The output file will contain a space-separated list of polymer concentrations _in units of Molar_ (see below), in the same order as in the CFE file. Please note that the concentrations can be in scientific notation like "4.47e-53" or "0.00e0", so you have to parse this properly.

Important: For systems of interest, the polymer basis can be quite large (hundreds of thousands of polymers). COFFEE can handle CFE files of this size. We need to make sure that our code can handle it as well.

### Optional temperature parameter
`tbnexplorer2` should take an optional `--temp <celsius>` parameter, which gets passed to `coffee-cli  --temp <celsius>`.

### Option to use Nupack 3 concentrations solver instead of COFFEE
Let's have an optional argument `--use-nupack-concentrations` to use Nupack's `concentrations` CLI tool instead of COFFEE to compute polymer concentrations. Path to `concentrations` is in .env as `NUPACK_CONCENTRATIONS_PATH`.

Command line: `concentrations -sort 0 -T {temperature_in_C} {base}`. The `-sort 0` flag preserves input polymer order between .ocx and .eq. The tool takes input files in a similar format as COFFEE called {base}.ocx (containing the matrix of polymers) and {base}.con (containing the concentrations of all the monomers in units of Molar). The .ocx file is a tab-delimited file with the first column being the polymer id (line number), the second column being always 1, and then the monomer composition of the polymer as for COFFEE, concluding in the free energy of that polymer.

After `concentration` is executed, {base}.eq is generated, which *adds* another column to the .ocx file (as last column) which contains the equilibrium concentrations of all the polymers in units of Molar.

## 4. Generate output TBN polymer matrix file
This file should be in the same format at the CFE file for `coffee-cli` described above, with an additional column for the concentration of that polymer computed by `coffee-cli`. 
The polymers should be **sorted** in order of decreasing concentration. 
There should also be some comments on top of the file indicating things like how many polymers there are in the polymer basis, and whether only partial computation was done (see below).
If the original input file was `<example>.tbn`, this file should be called `<example>.tbnpolymat`


# Partial computations
If the input .tbn file does not contain monomer concentrations, then we should do computation (1) above, and still output a .tbnpolymat file as in (4), except that it would be missing the column corresponding to concentrations. 
Additionally, we want to have command line options `--no-concentrations`, and `--no-free-energies` to disable the respective computations even if concentrations are given in the input .tbn file. 
Note that the `--no-free-energies` option also disables the concentrations computation since it's not possible to compute the concentrations without polymer free energies. Naturally, the `--no-free-energies` option also avoids adding the free energies column to the .tbnpolymat output file.


# Units
The concentration units are specified directly in the .tbn file using the `UNITS` keyword (see above). The monomer concentrations in the .tbn file are in the units specified by the UNITS line. These should be converted to Molar for COFFEE and then _back_ to the original units for the .tbnpolymat file. The comments on top of the .tbnpolymat file should also specify the units.


# User-friendly representation of polymers in .tbnpolys files
While .tbnpolymat files provide a matrix representation of polymers (and their concentrations), we also want a user-friendly (text-based) representation of polymers. This is provided by .tbnpolys files. We need to be able to output .tbnpolys files, as well as to take them as input (for functionality to be developed later). We should have a parses for these files in the Python package.

## Syntax of .tbnpolys files
- Comments are designated by "#" in the same way as for .tbn files.
- There should be at least one _empty line_ between different polymers (this is how we can tell when one polymer ends and the next begins). Note: Comments shouldn't count as empty lines
- Each polymer p is represented by its monomers, one per line
    - Each monomer can start with a multiplicity prefix of the form "n | " where n is the multiplicity of the monomer. If the multiplicity prefix is omitted, then we assume n = 1 by default.
    - After the multiplicity prefix, comes the monomer specification. This can either be the "name" of the monomer (as given in the .tbn file), or its binding site representation (eg: "a a* b2 b1"). Note: when the .tbnpolys file is given as user input, the binding sites can be given in any order, not necessarily in the same order as in the .tbn file. (Ie they should be mapped to vector representation of monomer.)


# Filtering output .tbnpolymat file and user-friendly presentation with `tbnexplorer2-filter`
The big picture is that often we want to know about which polymers certain monomers end up in at equilibrium. This can be hard to extract from the raw .tbnpolymat file. 
We make an additional command line tool `tbnexplorer2-filter` which takes a .tbn file as input. From the file name it infers the corresponding .tbnpolymat file as well. 

The `tbnexplorer2-filter` tool requires that the input .tbn file contains a `UNITS` keyword and monomer concentrations. If the .tbn file does not have `UNITS` specified, the tool will return an error, as it cannot function without concentration information.

The next thing on the command line is a space-separated list of _monomer names_: `<m1> <m2> ...`
The tool should output to the standard output only the polymers containing _all_ the monomers `m1 m2 ...`. If a monomer name repeats multiple times, we take this as the lower bound on the multiplicity of that monomer in the polymers to return. Order is ignored (we are looking at vectors only.)

If no monomer names are specified as input, then we do not do any filtering by monomer inclusion and return all polymers. (The output limits described below would still apply.)

## Output format
The output format should be user-friendly in `.tbnpolys` syntax. Polymer concentrations should be listed as comments after each polymer.
The order of the polymers should be in order of decreasing concentration.
The tool should also output what fraction (percent) of the _total concentration_ of all polymers in .tbnpolymat is the sum of the concentrations of the polymers matching the filtering criteria.
(To save space, we should not output the free energies of the complexes.)
Also, let's format the concentrations "nicely": i.e., instead of "9.99e+01 nM" we should say "99.9 nM"

## Output limit
Sometimes even with the filtering criteria, there are still too many polymers. 

First, we can limit the maximum number of polymers output with the optional `--num` parameter (short `-n`). This defaults to 100 if not explicitly specified.

Second, we add an optional command line argument `--min-concentration {conc}` where {conc} is a real non-negative number specifying the lower-bound on the concentration (in same units as \UNITS). 

Third, we add an optional command line argument `--percent-limit p` (short `-p`) where p is a percent number (real-value). 
The output should be restricted to those polymers whose concentration is above (p/100) fraction of the _total concentration_ of all polymers in the .tbnpolymat file. 

## Advanced filtering with constraints file
For more advanced filtering, the user can specify an optional `--constraints-file <filename>` argument. If this argument is given, no monomers should be specified on the command line and all the constraints should be only in the given constraints text file. 

The following constraints are allowed in the constraints file, one per line. As before, `<m1> <m2> ....` are monomer names. 
- `CONTAINS <m1> <m2> ...`: This is the same as specifying the monomers `<m1> <m2> ...` on the command line (see above).
- `EXACTLY <m1> <m2> ...`: This means only the polymer consisting exactly of the monomers `<m1> <m2> ...` and no other monomers. As with `CONTAINS`, multiplicity counts.
If multiple constraints are given (multiple lines), they should be treated as an OR. In other words, a polymer should be returned if it satisfies at least one of the constraints. The outputted polymers matching the constraints should always be in global order of decreasing concentration.


# Caching polymer basis
The most computationally intensive part of the pipeline is computing the polymer basis with Normaliz. Other parts, like using COFFEE, are typically much faster.
The polymer basis just depends on what the monomers are, not on their concentrations. Thus if we want to recompute polymer concentrations for new input monomer concentrations without changing what the monomers are, we should avoid re-computing the polymer basis. We do this as follows:

When `tbnexplorer2` is run with an input .tbn file, we create the A matrix and find the corresponding .tbnpolymat file (using the naming rules above) if it exists. If the .tbnpolymat has the keyword: \MATRIX-HASH: <hash>, we compare the hash of A with <hash>. If the hashes match, we can skip recomputing the polymer basis and instead load it from the .tbnpolymat file. Otherwise, we compute the polymer basis as normal and save the new hash. 

The keyword `\MATRIX-HASH: <hash>` should be somewhere close to the top of the .tbnpolymat file (without a comment marker prefix). The standard output of `tbnexplorer2` should indicate whether it had to re-generate the polymer basis or not (hashes matched).


# Advanced TBN specification using parametrized syntax (template syntax) for concentrations
The framework includes functionality for parsing .tbn files that include a parametrized specification for monomer concentrations.

## Command-line usage
To parse templated .tbn files, `tbnexplorer2` is passed an optional `--parametrized` argument followed by a list of numerical variable assignments such as `conc1=90.3 conc2=50`. 

## Template syntax
In the .tbn file, template variables and arithmetic expressions are specified using double curly braces: `{{expr}}`. These can be simple variable substitutions or arithmetic expressions involving variables. The expressions are evaluated using safe arithmetic evaluation and the results are used as concentration values.

### Supported operations
- Basic arithmetic: `+`, `-`, `*`, `/`
- Exponentiation: `**` (e.g., `{{2 ** 3}}` evaluates to 8)
- Parentheses for grouping: `{{(a + b) * c}}`
- Floating point numbers: `{{x * 1.5 + y * 0.5}}`

### Example .tbn file
```
\UNITS: nM
monomer1: a b*, {{conc1}}                     # Simple variable
monomer2: c d, {{conc2 + 10}}                 # Addition
monomer3: e f, 75.5                           # Literal value
monomer4: g h, {{base_conc * scale_factor}}   # Multiplication
monomer5: i j, {{(x + y) / 2}}                # Average calculation
```

## Parameter storage in .tbnpolymat
When a parametrized .tbn file is processed, the used parameter values are stored in the generated .tbnpolymat file using the `\PARAMETERS:` keyword in the header section. For example: `\PARAMETERS: conc1=90.3 conc2=50.0`

## Compatibility with tbnexplorer2-filter
The `tbnexplorer2-filter` command automatically reads parameters from the .tbnpolymat file when present. This ensures that filtering operations work correctly with parametrized .tbn files without requiring the user to re-specify the parameter values.


# Debugging Normaliz and 4ti2 tool calls
To aid debugging, we have the optional `--store-solver-inputs` parameter. This should create a `solver-inputs` subdirectory (if one doesn't exist) and copy all input files for Normaliz and 4ti2. This copy should occur prior to *each* call to Normaliz and 4ti2 (e.g., for computing the polymer basis, for computing irreducible canonical reactions in the extensions, etc.) The filenames should have the form `{base}...` where `{base}` is the input .tbn file name. The `...` should include information about which call this is (eg for computing the polymer basis or something else).
