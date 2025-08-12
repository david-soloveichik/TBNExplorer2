
# Background 
We are working with thermodynamics of abstract chemical systems (the model is called the TBN model). A TBN is defined in terms of binding sites (aka domains), monomers, and polymers (aka complexes). Conceptually, monomers are multisets of binding sites, and polymers are multisets of monomers. The big picture goal is given a set of monomers (described in terms of their binding sites) and their concentrations, we want to compute the equilibrium concentrations of the polymers. 

We will be building a Python library. The functionality should also be accessible via a command line (CLI) tool called `tbnexplorer2`.

# Input file encoding (.tbn) and TBN representation
A binding site is a symbol like "a" or "ahu34_8" or "5" (any string not containing white space nor any of the following: ",*|:" which will have special meaning). There are two types of binding sites, unstar and star. Star binding sites always have a "*" as the last character. The idea is that corresponding star and unstar binding sites bind (they are complementary).
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
monomer1: a a* b2 b1, 100
b2* b2* b1 b1  # this is a comment
a b2* b1, 50.7
```
Note that "#" followed by anything is a comment and should be ignored.
Monomers can have names which is indicated by a name following by ":" prior to the monomer specification.
Monomers can have concentrations indicated by a comma followed by the concentration after the monomer specification. There are two types of .tbn files: where each monomer has a concentration or none of the monomers have a concentration. 
Important: An error should be returned if some monomers have concentration and others not.

To avoid confusion, monomer names and binding sites should be distinct, and there should be an error otherwise. (Capitalization matters, so it's ok to have a binding site "c" and a monomer named "C", for example.) 
Monomer names should not have spaces or any of the prohibited symbols ",*|:"


# Polymers 
A polymer is a multiset of monomers. Note that we can have a polymer (or complex) consisting of just a single monomer (we would call this a "singleton polymer").

Internally, we will represent a polymer as vector of non-negative integers corresponding to the counts of the monomers in it. Given our representation of monomers in matrix A above, if x is a polymer, then b = A.x is a vector of binding site *excesses* in the polymer x. If b = 0, then every binding site is paired to its complement.


**Star-limiting restriction**: We always enforce the "star-limited restriction" described below, and give an error if it is not satisfied. Intuitively, a TBN is star-limited if for every binding site there is at least as much x as x* (totalled over all monomers). More formally:

- When monomer concentrations are given in the TBN file: Let c be the vector of the concentrations of the monomers. Then the TBN is *star-limited* if A.c ≥ 0 (component-wise).
- When monomer concentrations are not given, then c = Transpose[(1,1,...,1)] for the star-limited restriction check.


# Additional Utilities
- We want to have Python utility functions for running `Normaliz` to solve Hilbert basis problems and retrieving its output. It is a bit tricky to extract the Hilbert basis from the command line output of Normaliz. Important: please see how `/Users/dsolov/Documents/Development/VibeDevelopment/TBNCanonicalReactionsEnumerator` retrieves the Hilbert basis from the output of Normaliz.
The path to Normaliz should be hardcoded in the code, but easily changeable (as a constant on top of the relevant file).


# Computations over TBNs
Given an input TBN, we will be interested in computing several things. 

## 1. Polymer basis (aka Hilbert basis of polymers)
The polymer basis is conceptually the set of "unsplittable" polymers, the polymers that cannot be decomposed into two without losing some bonding. We compute it by solving a special linear problem. 

First, we need to make sure that for every binding site x, there is a column in A consisting of all zeros and a "-1" in the position corresponding to x. This is intuitively equivalent to having a singleton monomer {x*} in the TBN. If there are some binding sites for which such a singleton monomer does not exist in the TBN, then we add these columns to A as if we had more monomers. Importantly, we need to keep track of how many monomers we had originally, so that we can remove these fake added monomers later. Call this modified matrix A'. 

Now we want to find the Hilbert basis of the cone corresponding to solutions of A'.x = 0. We do this using `Normaliz`. Let H be the Hilbert basis output by Normaliz for this problem. To get the *polymer basis*, we remove the last entries of vectors in H corresponding to the "fake" singleton monomers (that we added to form A'), and remove any resulting duplicates in H. 

Important: The Hilbert basis H may be very large! There could be hundreds of thousands of vectors in H. Thus these operations need to be very efficient. 

We optionally save the polymer basis into a user-friendly text file. If the command line option `--user-friendly-polymer-basis` is given to out CLI tool `tbnexplorer2`, it should do the following: 
If given [example].tbn as input, it should create a file called [example]-polymer-basis.txt. This file should represent the polymers in the polymer basis in the following user-friendly way:
- There should be an empty line between polymers
- Each polymer p is represented by its monomers, one per line
    - The line should start with "n | " where n is the multiplicity of the monomer in p
    - After this, the monomer is represented as just its "name" if this name was given in the input .tbn file.
    - If the monomer's name was not given, then the monomer should be represented as its binding site representation; eg: "a a* b2 b1"
    Note: the binding sites should appear in the same order as in the original .tbn input file (but extra whitespace should be removed).

The CLI tool should return to the command line the number of polymers in the polymer basis.

## 2. Compute polymer free energies
Given the polymers in the polymer basis, we want to compute their "free energies". 
In general, we want to implement the following: Given a matrix of polymers, with a polymer per row (ie vector of the monomer counts of that polymer), we want to compute the "free energies" of each polymer. 

Intuitively, the free energy of a polymer x is dG(x) = -[number of bonds in that polymer]. 
We exclude self-binding within a monomer from that calculation (e.g., in a monomer like {a, b, b*}, we don't count the b-b* bond). This is because this bond is always there.

We can compute the number of bonds in polymer x as: (Sum[|A|.x] - Sum[A.x])/2, where |A| is the same as A but with absolute value applied to all entries, and Sum[v] sums all the elements of v. 
Thus, intuitively, Sum[|A|.x] is the total number of binding sites in x (excluding self-binding within a monomer), and Sum[A.x] is the total excess of unstar binding sites. Subtracting the two gives twice the number of bonds formed (since each bond involves exactly 2 binding sites), so we divide by 2. 

## 3. Compute equilibrium polymer concentrations
The big picture is that we want to compute the equilibrium concentrations of all the polymers in the polymer basis. We use the command line tool COFFEE for this, with the executable `/Users/dsolov/Documents/ResearchTools/coffee/crates/coffee-cli/target/release/coffee-cli`.

Please see `/Users/dsolov/Documents/ResearchTools/coffee/README.md` for COFFEE documentation. Roughly, it takes two input files: CFE and CON. 
- The CFE file contains the following matrix: Each line corresponds to a polymer in the polymer basis (i.e., counts of monomers in that polymer), followed by the free energy of that polymer. 
- The CON file contains the concentration of each monomer, one per line. So there should be as many lines as columns in the CFE file minus one (the free energy). The monomer concentrations should be as given in the input .tbn file, _in units of Molar_ (see below).
Use the `-o [filename]` flag for `coffee-cli` to write its output in [filename]. The output file will contain a space-separated list of polymer concentrations _in units of Molar_ (see below), in the same order as in the CFE file. Please note that the concentrations can be in scientific notation like "4.47e-53" or "0.00e0", so you have to parse this properly.

Important: For systems of interest, the polymer basis can be quite large (hundreds of thousands of polymers). COFFEE can handle CFE files of this size. We need to make sure that our code can handle it as well.

## 4. Generate output TBN polymer matrix file
This file should be in the same format at the CFE file for `coffee-cli` described above, with an additional column for the concentration of that polymer computed by `coffee-cli`. 
The polymers should be **sorted** in order of decreasing concentration. 
There should also be some comments on top of the file indicating things like how many polymers there are in the polymer basis, and whether only partial computation was done (see below).
If the original input file was `[example].tbn`, this file should be called `[example].tbnpolymat`


# Partial computations
If the input .tbn file does not contain monomer concentrations, then we should do computation (1) above, and still output a .tbnpolymat file as in (4), except that it would be missing the column corresponding to concentrations. 
Additionally, we want to have command line options `--no-concentrations`, and `--no-free-energies` to disable the respective computations even if concentrations are given in the input .tbn file. 
Note that the `--no-free-energies` option also disables the concentrations computation since it's not possible to compute the concentrations without polymer free energies. Naturally, the `--no-free-energies` option also avoids adding the free energies column to the .tbnpolymat output file.


# Units
Our `tbnexplorer2` tool should take an additional command line argument `--concentration-units` which should default to "nM" (nanoMolar) if not given. Other options are "pM (picoMolar), "uM" (microMolar), "mM" (milliMolar), and "M" (Molar).
The input .tbn file specifies the monomer concentrations in the units given by the `--concentration-units` parameter. This should be converted to Molar for COFFEE and then _back_ to the desired units for the .tbnpolymat file. The comments on top of the .tbnpolymat file should also specify the units.


# Filtering output .tbnpolymat file
The big picture is that often we want to know about which polymers certain monomers end up in at equilibrium. This can be hard to extract from the raw .tbnpolymat file.

We make an additional command line tool `tbnexplorer2-filter` which takes a .tbn file as input. From the file name it infers the corresponding .tbnpolymat file as well. 

The next thing on the command line is a space-separated list of _monomer names_: `m1 m2 ...`

The tool should output to the standard output only the polymers containing _all_ the monomers `m1 m2 ...`. If a monomer name repeats multiple times, we take this as the lower bound on the multiplicity of that monomer in the polymers to return.

If no monomer names are specified as input, then we do not do any filtering by monomer inclusion and return all polymers. (The output limits described below would still apply.)

## Output format
The output format should be user-friendly like in the `[example]-polymer-basis.txt` file, except that polymer concentrations are listed. 
The order of the polymers should be in order of decreasing concentration.
The tool should also output what fraction (percent) of the _total concentration_ of all polymers in .tbnpolymat is the sum of the concentrations of the polymers matching the filtering criteria.
To save space, we should not output the free energies of the complexes.
Also, let's format the concentrations "nicely": i.e., instead of "9.99e+01 nM" we should say "99.9 nM"

## Output limit
Sometimes even with the filtering criteria, there are still too many polymers. 

First, we can limit the maximum number of polymers output with the optional `--num` parameter (short `-n`). This defaults to 100 if not explicitly specified.

Second, we add an optional command line argument `--percent-limit p` (short `-p`) where p is a percent number (real-value). 
The output should be restricted to those polymers whose concentration is above (p/100) fraction of the _total concentration_ of all polymers in the .tbnpolymat file. 
