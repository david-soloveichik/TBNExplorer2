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

We use Normaliz (or 4ti2 with `--use-4ti2`) to get the Hilbert basis of this cone. We call these the **irreducible canonical reactions**. 


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
- Repeat until no irreducible canonical reactions remain. 
End iterations
After finishing the iterations, any remaining polymers p with `μ(p) = 0` (i.e., whose concentration exponents are still not assigned) are removed from consideration. (Intuitively, they are the off-target polymers that cannot be generated by any canonical reaction from on-target polymers and we don't care about them.)

Important: this algorithm must be implemented very efficiently since there might be hundreds of thousands of irreducible canonical reactions and hundreds of thousands of polymers.

## Output text file showing irreducible canonical reactions
If an optional command line argument `--output-canonical-reactions` is given:  
Output a .txt file showing all irreducible canonical reactions. They should be ordered by the iteration of the IBOT algorithm in which they were part of R. The μ_min should be recorded. Also the polymers that get assigned `μ(p) = μ_min` should be indicated somehow in the reaction in a compact way.

## Output .tbnpolys file
Generate a .tbnpolys file containing all polymers in the polymer basis (i.e., on-target and off-target). For each, indicate its concentration exponent as `μ: {value}` on a separate line following the polymer.

The on-target polymers should come first, and the off-target afterward. There should be comments indicating the on-target group and the off-target group. The off-target polymers should be *ordered* by their concentration exponent in increasing order (i.e., smallest to largest). 

## Output .tbn file
If given the optional argument `--generate-tbn {c} {units}`, we make a .tbn file as follows.
- It has \UNITS as specified by {units}
- It has the same monomers as the original input .tbn file
- Each monomer i is assigned concentration `... + p[i] * ((c'/ρH20) ** μ(p))*ρH20 + ...` where:
    - `c'` is `c` converted to Molar from {units}
    - `ρH20 = 55.14` is the density of water at 37C. (Intuitively, the process of dividing and them multiplying by `ρH20` is because the IBOT algorithm logic works in units of "mole fractions" rather than Molar.)
    - the sum is over all polymers p where monomer i appears with (non-zero) multiplicity `p[i]`.
    - `μ(p)` is the concentration exponent of polymer p computed in the IBOT algorithm.


# Testing the big picture

Please generate a separate python file which tests the following.

To close the loop and test the pipeline:
Suppose we start with some `input.tbn` file (without concentrations), and execute `tbnexplorer2-ibot` to get the concentration exponents `μ(p)` for all polymers (in the polymer basis). 
We also use `--generate-tbn {c} nM` to generate `output.tbn` file (with concentrations).
Next run `tbnexplorer2 output.tbn` (and possibly `tbnexplorer2-filter output.tbn` for nicer output) to get polymer concentrations from `output.tbn` (again for all polymers in the polymer basis). 
Then the concentration of each polymer p thus obtained should be the same (or similar due to floating point errors) to `((c'/ρH20) ** μ(p))*ρH20`, where `μ(p)` was the concentration coefficient generated by `tbnexplorer2-ibot`.

If the results are quantitatively off, there might be deep conceptual issues that need to be addressed separately, but it would be good to know how far off they are.


# Bounding specific off-target polymer concentrations (undesired polymers)

Often we are only interested in an upper bound on the concentrations of some specific set undesired off-target polymers p_1,...,p_k.
In particular, for large systems, the generation of all irreducible canonical reactions (by Normaliz or 4ti2) takes too long or is completely infeasible.
Rather than generating the exact equilibrium concentrations via the IBOT algorithm, we can much more more efficiently compute an upper bound on the p_i by narrowing our focus to reactions that directly produce some p_i instead of examining the full set of irreducible canonical reactions.  

## Mathematical background

Let P = all canonical reactions producing a given undesired polymer p.
We want to generate a subset T of P with the following properties:

1. Every x \in P can be written as:
    u + v
    where u \in T,
    and v is some _canonical reaction_ (possibly zero)
2. T is minimal in the sense that:
    For every u \in T, there do not exist a canonical reaction u' that produces p
    and some canonical reaction v \neq 0 such that:
    u = u' + v.

The way we find such a T for each undesired polymer p_i (separately) is via the following linear problem:

For each p_i we do the following:
Construct the linear system:
B * r = 0
S * r ≥ 0
e_i * r > 0  (equivalently, e_i * r ≥ 1 over integers)
where B and S are as described above, and e_i selects p_i in r. 
We use 4ti2's `zsolve` to compute the minimal inhomogeneous solutions (module generators) for this system with strict inequality. Call this T_i.
Repeat for next i.

I claim that each T_i satisfies the two properties above.

Our reduced set of canonical reactions is now the union of the T_i. When we compute `μ(pi)` using this reduced set of reactions, we will get something that is a lower-bound on the true `μ(pi)`. This gives us an upper-bound on `(c'/ρH20)^μ(pi)` since mole fraction `c'/ρH20` is always less than 1.

## Implementation Note: Variable Splitting Instead of Explicit S Matrix

The implementation uses a **variable splitting technique** rather than explicitly constructing the S matrix described above. This mathematically equivalent approach:

1. **Splits on-target polymer variables** into positive and negative parts (r = r_pos - r_neg), allowing them to take any sign
2. **Keeps off-target polymer variables** as single non-negative variables (r ≥ 0)

This transformation implicitly enforces S*r ≥ 0 (no off-target reactants) because:

- Off-target polymers are constrained to r ≥ 0, so they can only appear as products (positive values)
- On-target polymers can be positive or negative, allowing them to be reactants or products

**Why this approach?** Variable splitting is more efficient as it:

- Avoids constructing the large S matrix explicitly
- Reduces the number of inequality constraints
- Is a standard technique in polyhedral computation that solvers like Normaliz handle efficiently

To implement this functionality, we introduce an optional command line argument `--upper-bound-on-polymers {undesired_off_target}.tbnpolys` where the file is in .tbnpolys syntax specifying the p_i.

When this option is used, we cannot use the option `--generate-tbn {c} {units}` since we do not know all the `μ(p)` for all off-target p.
The upper bounds computation itself always uses 4ti2's `zsolve` (regardless of the `--use-4ti2` flag), while `--use-4ti2` controls whether to use 4ti2 or Normaliz for the polymer basis computation.
We can still use `--output-canonical-reactions` to show the relevant information for the reactions we generated in the new reduced way.

To test the bounding functionality, we can add _all_ off-target polymers (that can be produced by some reaction from on-target polymers) to `{undesired_off_target}.tbnpolys`. Then if our implementation is correct, the concentration coefficients obtained should be _identical to_ (not just lower-bound) the case without `--upper-bound-on-polymers` option. As the test system for this, we can use files in `extensions/my_inputs/testing_bounding_IBOT/`:

- TBN system: `and_gate_noA.tbn`
- On-target polymers: `and_gate_noA_on-target.tbnpolys`
- All off-target polymers (that can be produced by some reaction from on-target polymers): `and_gate_noA_all-off-target.tbnpolys`
