Extensions are not part of the main functionality of TBNExplorer2, but rely on the framework to perform related computation. The code for the extensions is in the `extensions/` subfolder for separation of concerns. 

# Extension: Enumerating Irreducible Canonical Reactions

## On-target and off-target polymers
The big picture is that we are given a .tbn file (without concentrations, i.e., no \UNITS) from which we compute the polymer basis as usual. This is the set of _all_ polymers.
We are also given an input .tbnpolys file listing the so-called *on-target* polymers. The rest of the polymers are considered off-target.
There should be an error if the input .tbnpolys file contains a polymer that is not in the polymer basis.

## Mathematical Background
Suppose there are n polymers and m monomers. 

### Polymer reactions
We consider reactions *between polymers*. A reaction r corresponds to a vector of integers of length n: 
r[p] = the change of count of polymer p in reaction r 
Examples:
1. If the reactants of r include p with multiplicity 1, and the products do not include p, then r[p] = -1.
2. If the reactants of r include p with multiplicity 1, and the products include p with multiplicity 3, then r[p] = 2.
3. If the reactants of r no not include p, and the products include p with multiplicity 2, then r[p] = 2.

Important note: We will actually not allow "catalysts" (i.e. example 2 where p appears both as a reactant and product), so our representation as a vector can be _uniquely_ decoded in terms of reactants and products. 

### B matrix and S matrix; valid canonical reactions
Important note: We will be dealing with potentially hundreds of thousands of polymers so we need to have efficient data structures.

Matrix B (m×n): B[i,p] = count of monomer i in polymer p. 
Any *valid* reaction r must satisfy the conservation of mass of monomers: `B*r = 0`.

Matrix S ((n - num_on_target) x n): Selects the off-target polymers. Specifically, `S*r` is the vector of off-target polymer components of r.

We say a reaction r is *canonical* if it is valid and `S*r ≥ 0`. In other words, no off-target polymer appears as a reactant. Intuitively, canonical reactions are those that create off-target polymers from purely on-target ones, and these are the reactions we must be worried about. 

### Linear problem, Hilbert basis, and irreducible canonical reactions
Based on the above, the space of all canonical reactions corresponds to the integer solutions (r) of the following system, which is a "cone":
B*r = 0
S*r ≥ 0

Now we are interested in the irreducible reactions in the space of all canonical reactions: that cannot be written as a sum of two other canonical reactions. These correspond exactly to the Hilbert basis of the above cone. 

We use Normaliz (or 4ti2) to get the Hilbert basis of this cone. We call these the **irreducible canonical reactions**. 


# Extension: Iterative Balancing of Off-Target Polymers ("IBOT")
This extension implements an iterative algorithm for assigning concentrations to off-target polymers (such that they are in detailed balance). It depends on the enumeration of irreducible canonical reactions by the previous extension.

The functionality of this extension can be accessed with the command line tool is called `tbnexplorer2-ibot`.

## Checking detailed balance of on-target polymers
We make sure that all irreducible canonical reactions that are entirely over on-target polymers (i.e., don't have any off-target polymers as products) have the same number of reactants as products (including their multiplicity). For example, if all the Pi's are on-target, then this is ok:
P1 + 2 P2 -> P3 + P4 + P5
This is not ok because there are more products than reactants:
P1 + P2 -> P3 + 2 P5
If there is a reaction which violates this, we return an error: "On-target polymers not in detailed balance" and show the violating reaction. 

## Computing imbalance and novelty of irreducible canonical reactions
The "IBOT" algorithm described below iteratively assigns "concentration exponents" to off-target polymers. For a polymer p, its concentration exponent is written as `μ(p)`. 

All on-target polymers p have concentration exponents pre-assigned at `μ(p) = 1`.
Prior to being assigned their concentration exponent, an off-target polymer p has `μ(p) = 0`.

Suppose we have already assigned concentration exponents to some set of off-target polymers. Given an irreducible canonical reaction r, we are interested in two key quantities:

1. novelty `l(r)`: how many off-target polymers appear in r that we haven't assigned concentration exponents to yet

2. imbalance `k(r)`: Let μ1 be the sum of the concentration exponents of the reactants of r weighed by their multiplicity. Let μ2 be the sum of the concentration exponents of the products of r weighed by their multiplicity. Then the imbalance of r is `k(r) = μ1 - μ2`.

Note: We can compute μ1 as the dot product of the negative coordinates p of r with `-μ(p)`. Similarly, we can compute μ2 as the dot product of the positive coordinates p of r with `μ(p)`.

## IBOT algorithm
In each iteration:
- We eliminate all irreducible canonical reactions r with `l(r) = 0`.
- For each remaining irreducible canonical reaction r we compute the imbalance-novelty ratio `k(r)/l(r)`.
- We find the reactions with the smallest imbalance-novelty ratio. Call this smallest imbalance-novelty ratio μ_min and call these reactions R. 
- All polymers p appearing in any reaction r in R such that `μ(p) = 0` (i.e., those that have not previously been assigned a concentration exponent) get assigned `μ(p) = μ_min`
- Repeat until no irreducible canonical reactions remain. This is the same as saying that all off-target polymers have been assigned concentration exponents.

Important: this algorithm must be implemented very efficiently since there might be hundreds of thousands of irreducible canonical reactions and hundreds of thousands of polymers.

## Output .tbnpolys file
Generate a .tbnpolys file containing all polymers in the polymer basis (i.e., on-target and off-target). For each, indicate its concentration exponent as `μ: {value}` on a separate line following the polymer.

The on-target polymers should come first, and the off-target afterward. There should be comments indicating the on-target group and the off-target group. The off-target polymers should be *ordered* by their concentration exponent in increasing order (i.e., smallest to largest). 

## Output .tbn file
If given the optional argument `--generate-tbn {c} {units}`, we make a .tbn file as follows.
- It has \UNITS as specified by {units}
- It has the same monomers as the original input .tbn file
- Each monomer i is assigned concentration `... + p[i] * c ** μ(p) + ...` where:
    - the sum is over all polymers p where monomer i appears with (non-zero) multiplicity `p[i]`.
    - `μ(p)` is the concentration exponent of polymer p computed in the IBOT algorithm.


# Testing the big picture
Please generate a separate python file which tests the following.

To close the loop and test the pipeline:
Suppose we start with some `input.tbn` file (without concentrations), and execute `tbnexplorer2-ibot` to get the concentration exponents `μ(p)` for all polymers (in the polymer basis). 
We also use `--generate-tbn {c} nM` to generate `output.tbn` file (with concentrations).
Next run `tbnexplorer2 output.tbn` (and possibly `tbnexplorer2-filter output.tbn` for nicer output) to get polymer concentrations from `output.tbn` (again for all polymers in the polymer basis). 
Then the concentration of each polymer p thus obtained should be the same (or similar due to floating point errors) to `c ** μ(p)`, where `μ(p)` was the concentration coefficient generated by `tbnexplorer2-ibot`.

If the results are quantitatively off, there might be deep conceptual issues that need to be addressed separately, but it would be good to know how far off they are.