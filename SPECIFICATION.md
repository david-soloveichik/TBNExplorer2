
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

We need to be able to save the polymer basis into a text file. The CLI tool `tbnexplorer2`, if given [example].tbn as input, should create a file called [example]-polymer-basis.txt. This file should represent the polymers in the polymer basis in the following user-friendly way:
- There should be an empty line between polymers
- Each polymer p is represented by its monomers, one per line
    - The line should start with "n | " where n is the multiplicity of the monomer in p
    - After this, the monomer is represented as just its "name" if this name was given in the input .tbn file.
    - If the monomer's name was not given, then the monomer should be represented as its binding site representation; eg: "a a* b2 b1"
    Note: the binding sites should appear in the same order as in the original .tbn input file (but extra whitespace should be removed).

The CLI tool should return to the command line the number of polymers in the polymer basis.

## 2. FOR NOW WE END WITH (1), TO BE CONTINUED
In particular, we'll deal with monomer concentrations further in later functionality.

