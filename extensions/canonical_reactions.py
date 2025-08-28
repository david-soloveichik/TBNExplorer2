"""Module for computing irreducible canonical reactions.

This module provides functionality to enumerate irreducible canonical reactions
between on-target and off-target polymers in a TBN system.
"""

from pathlib import Path
from typing import List, Optional, Set, Tuple

import numpy as np

from tbnexplorer2.fourtitwo import FourTiTwoRunner
from tbnexplorer2.model import TBN
from tbnexplorer2.normaliz import NormalizRunner
from tbnexplorer2.tbnpolys_io import TbnpolysParser


class Reaction:
    """Represents a reaction between polymers."""

    def __init__(self, vector: np.ndarray, polymer_names: Optional[List[str]] = None):
        """
        Initialize a reaction.

        Args:
            vector: Reaction vector where negative values are reactants, positive are products
            polymer_names: Optional list of polymer names for display
        """
        self.vector = vector
        self.polymer_names = polymer_names

    def get_reactants_and_products(self) -> Tuple[List[Tuple[int, int]], List[Tuple[int, int]]]:
        """
        Get reactants and products from reaction vector.

        Returns:
            Tuple of (reactants, products) where each is a list of (polymer_index, multiplicity)
        """
        reactants = []
        products = []

        for i, count in enumerate(self.vector):
            if count < 0:
                reactants.append((i, abs(count)))
            elif count > 0:
                products.append((i, count))

        return reactants, products

    def is_balanced(self) -> bool:
        """
        Check if reaction has same number of reactants and products (including multiplicity).

        Returns:
            True if sum of reactant multiplicities equals sum of product multiplicities
        """
        reactants, products = self.get_reactants_and_products()
        reactant_count = sum(mult for _, mult in reactants)
        product_count = sum(mult for _, mult in products)
        return reactant_count == product_count

    def __str__(self) -> str:
        """String representation of the reaction."""
        reactants, products = self.get_reactants_and_products()

        def format_side(polymers):
            terms = []
            for idx, mult in polymers:
                name = self.polymer_names[idx] if self.polymer_names else f"P{idx}"
                if mult == 1:
                    terms.append(name)
                else:
                    terms.append(f"{mult} {name}")
            return " + ".join(terms) if terms else "0"

        return f"{format_side(reactants)} -> {format_side(products)}"


class CanonicalReactionsComputer:
    """Computes irreducible canonical reactions for a TBN system.

    A canonical reaction is one where off-target polymers can only appear as products,
    not as reactants. Mathematically, this is expressed as S*r >= 0, where S is a
    selection matrix for off-target polymers and r is the reaction vector.

    IMPORTANT IMPLEMENTATION NOTE:
    Rather than explicitly constructing and using the S matrix, this implementation
    enforces the canonical constraint through variable splitting:
    - On-target polymers are split into positive and negative variables (can be +/- in reactions)
    - Off-target polymers use single non-negative variables (can only be >= 0 in reactions)
    This approach implicitly enforces S*r >= 0 without needing the S matrix explicitly.
    """

    def __init__(self, tbn: TBN, use_4ti2: bool = False):
        """
        Initialize the canonical reactions computer.

        Args:
            tbn: The TBN model
            use_4ti2: Whether to use 4ti2 instead of Normaliz
        """
        self.tbn = tbn
        self.use_4ti2 = use_4ti2
        self.polymers = None
        self.on_target_indices = None
        self.off_target_indices = None
        self.B_matrix = None
        self.S_matrix = None

    def load_on_target_polymers(self, tbnpolys_file: Path, polymer_basis: List[np.ndarray]) -> Set[int]:
        """
        Load on-target polymers from .tbnpolys file and identify their indices.

        Args:
            tbnpolys_file: Path to .tbnpolys file containing on-target polymers
            polymer_basis: List of all polymers in the basis (as monomer count vectors)

        Returns:
            Set of indices of on-target polymers in the polymer basis

        Raises:
            ValueError: If an on-target polymer is not in the polymer basis
        """
        # Parse on-target polymers
        parser = TbnpolysParser(self.tbn)
        on_target_polymers_raw = parser.parse_file(tbnpolys_file)

        # Convert to monomer count vectors
        on_target_polymers = []
        for polymer_raw in on_target_polymers_raw:
            # Create monomer count vector
            counts = np.zeros(len(self.tbn.monomers), dtype=int)
            for multiplicity, monomer in polymer_raw:
                # Find monomer index
                monomer_idx = self.tbn.monomers.index(monomer)
                counts[monomer_idx] += multiplicity
            on_target_polymers.append(counts)

        # Find indices of on-target polymers in the polymer basis
        on_target_indices = set()
        for on_target in on_target_polymers:
            found = False
            for i, polymer in enumerate(polymer_basis):
                if np.array_equal(on_target, polymer):
                    on_target_indices.add(i)
                    found = True
                    break

            if not found:
                raise ValueError(f"On-target polymer {on_target} not found in polymer basis")

        return on_target_indices

    def setup_matrices(self, polymer_basis: List[np.ndarray], on_target_indices: Set[int]):
        """
        Set up B and S matrices for the canonical reactions computation.

        Args:
            polymer_basis: List of all polymers (as monomer count vectors)
            on_target_indices: Set of indices of on-target polymers
        """
        self.polymers = polymer_basis
        self.on_target_indices = on_target_indices
        self.off_target_indices = set(range(len(polymer_basis))) - on_target_indices

        n_monomers = len(self.tbn.monomers)
        n_polymers = len(polymer_basis)
        n_off_target = len(self.off_target_indices)

        # B matrix: B[i,p] = count of monomer i in polymer p
        # Shape: (n_monomers, n_polymers)
        self.B_matrix = np.zeros((n_monomers, n_polymers), dtype=int)
        for p, polymer in enumerate(polymer_basis):
            self.B_matrix[:, p] = polymer

        # S matrix: Selects off-target polymers
        # Shape: (n_off_target, n_polymers)
        self.S_matrix = np.zeros((n_off_target, n_polymers), dtype=int)
        for i, p in enumerate(sorted(self.off_target_indices)):
            self.S_matrix[i, p] = 1

    def compute_irreducible_canonical_reactions(self) -> List[Reaction]:
        """
        Compute the irreducible canonical reactions.

        These are reactions that:
        1. Conserve mass (B*r = 0)
        2. Are canonical (S*r >= 0, no off-target reactants)
        3. Are irreducible (cannot be written as sum of two other canonical reactions)

        IMPORTANT: The canonical constraint S*r >= 0 is enforced implicitly through
        variable splitting rather than constructing an explicit S matrix. See detailed
        comments in the implementation below.

        Returns:
            List of Reaction objects representing irreducible canonical reactions
        """
        if self.B_matrix is None or self.S_matrix is None:
            raise RuntimeError("Matrices not set up. Call setup_matrices first.")

        n_polymers = self.B_matrix.shape[1]

        # Construct the cone system for Normaliz/4ti2:
        # We want: B*r = 0 and S*r >= 0
        # This is equivalent to finding the Hilbert basis of the cone:
        # { r : B*r = 0, S*r >= 0 }

        # For Normaliz/4ti2, we need to express this as:
        # Equations: B*r = 0
        # Inequalities: S*r >= 0

        # However, Normaliz expects the format: A*x = 0, x >= 0
        # So we need to transform our system

        # We'll use a standard transformation:
        # Let r = r+ - r- where r+, r- >= 0
        # Then B*r = B*(r+ - r-) = 0 becomes B*r+ - B*r- = 0
        # And S*r = S*(r+ - r-) >= 0 becomes S*r+ - S*r- >= 0

        # But for canonical reactions, we want a more specific structure.
        # We'll compute the Hilbert basis of the system directly.

        # Create augmented system for equations and inequalities
        # We need to handle this carefully for large systems

        # Combine equations and inequalities into a single system for Normaliz
        # The system is: find r such that B*r = 0 and S*r >= 0

        # Use the dual cone formulation
        # We create a matrix that encodes both constraints

        # For efficiency with large systems, we'll use a specialized approach
        # Create the constraint matrix for the cone
        # A_eq = self.B_matrix  # Equations: B*r = 0 (kept for reference)
        # A_ineq = -self.S_matrix  # Inequalities: -S*r <= 0 (i.e., S*r >= 0) (kept for reference)

        # Combine into single system for Normaliz
        # We'll use the inhomogeneous system approach

        # Create runner
        runner = FourTiTwoRunner() if self.use_4ti2 else NormalizRunner()

        # For the cone { r : B*r = 0, S*r >= 0 }, we need to be careful
        # We'll create an augmented system that Normaliz can handle

        # Create the full constraint matrix
        # Stack B and -S vertically

        # For Normaliz, we need to encode this as a cone problem
        # We'll use slack variables for inequalities

        # Augmented system: [B 0; S I] * [r; s] = [0; 0], s >= 0
        # This gives us B*r = 0 and S*r + s = 0, with s >= 0
        # Which means S*r = -s <= 0, so -S*r >= 0, so S*r >= 0 (after sign flip)

        # Actually, for canonical reactions, we want a simpler approach
        # We'll compute the Hilbert basis of the intersection of:
        # 1. The kernel of B (mass conservation)
        # 2. The non-negative orthant for off-target components

        # This is best done by computing the Hilbert basis of the lifted cone

        # Create the matrix for the lifted cone problem
        # Variables: r (reactions), one for each polymer
        # Constraints: B*r = 0 (mass conservation)
        #             r[i] >= 0 for i in off_target_indices (canonicality)

        # We need to handle the fact that on-target polymers can have negative values
        # So we'll use a different encoding

        # VARIABLE SPLITTING TO ENFORCE CANONICAL CONSTRAINT:
        # ====================================================
        # Instead of explicitly using S*r >= 0, we enforce the canonical constraint
        # (no off-target reactants) through variable transformation:
        #
        # Split variables into positive and negative parts for on-target polymers only
        n_on_target = len(self.on_target_indices)
        n_off_target = len(self.off_target_indices)

        # Variable order in lifted space: [r_on_target_pos, r_on_target_neg, r_off_target]
        # Total variables: 2 * n_on_target + n_off_target
        #
        # Key insight: By keeping off-target variables as single non-negative variables,
        # we ensure they can never be negative (reactants) in the original reaction vector.

        # Create the new B matrix for the split variables
        on_target_list = sorted(self.on_target_indices)
        off_target_list = sorted(self.off_target_indices)

        B_lifted = np.zeros((self.B_matrix.shape[0], 2 * n_on_target + n_off_target), dtype=int)

        # Fill in the lifted B matrix
        for i, p in enumerate(on_target_list):
            # On-target: split into positive and negative parts
            # Original: r[p] = pos - neg, allowing any sign
            B_lifted[:, i] = self.B_matrix[:, p]  # Positive part
            B_lifted[:, n_on_target + i] = -self.B_matrix[:, p]  # Negative part (subtract)

        for i, p in enumerate(off_target_list):
            # Off-target: single non-negative variable
            # Original: r[p] = h[i] >= 0, ensuring products only
            B_lifted[:, 2 * n_on_target + i] = self.B_matrix[:, p]

        # Compute Hilbert basis of { x >= 0 : B_lifted * x = 0 }
        hilbert_basis = runner.compute_hilbert_basis(B_lifted)

        if not hilbert_basis:
            return []

        # Convert back to reaction vectors
        reactions = []
        for h_vector in hilbert_basis:
            # Reconstruct the original reaction vector from lifted space
            reaction = np.zeros(n_polymers, dtype=int)

            # On-target polymers: r[p] = h_pos[i] - h_neg[i]
            # Can be positive (product), negative (reactant), or zero
            for i, p in enumerate(on_target_list):
                reaction[p] = h_vector[i] - h_vector[n_on_target + i]

            # Off-target polymers: r[p] = h[i] >= 0
            # Always non-negative, ensuring they only appear as products
            # This is how we enforce S*r >= 0 without explicit S matrix
            for i, p in enumerate(off_target_list):
                reaction[p] = h_vector[2 * n_on_target + i]

            # Skip trivial reactions (all zeros)
            if np.any(reaction != 0):
                reactions.append(Reaction(reaction))

        return reactions

    def compute_irreducible_canonical_reactions_for_targets(self, target_polymer_indices: Set[int]) -> List[Reaction]:
        """
        Compute irreducible canonical reactions that produce specific target polymers.

        This is used for computing upper bounds on specific off-target polymer concentrations.
        For each target polymer p_i, we compute module generators for the system:
        - B*r = 0 (mass conservation)
        - S*r >= 0 (canonical reactions - no off-target reactants)
        - e_i*r >= 1 (must produce target polymer p_i)

        IMPORTANT: The S*r >= 0 constraint is enforced implicitly through variable splitting:
        - On-target polymers are split into positive and negative parts (can be reactants or products)
        - Off-target polymers are kept as single non-negative variables (can only be products)
        This ensures off-target polymers never appear as reactants, which is what S*r >= 0 enforces.

        Args:
            target_polymer_indices: Set of polymer indices to target (must be off-target)

        Returns:
            List of Reaction objects that produce at least one target polymer

        Raises:
            ValueError: If target polymers are not all off-target
            RuntimeError: If matrices not set up
        """
        if self.B_matrix is None or self.S_matrix is None:
            raise RuntimeError("Matrices not set up. Call setup_matrices first.")

        # Validate that all target polymers are off-target
        invalid_targets = target_polymer_indices & self.on_target_indices
        if invalid_targets:
            raise ValueError(f"Target polymers must be off-target. Invalid indices: {invalid_targets}")

        # Check that target polymers exist
        n_polymers = self.B_matrix.shape[1]
        invalid_indices = {idx for idx in target_polymer_indices if idx >= n_polymers}
        if invalid_indices:
            raise ValueError(f"Target polymer indices out of range: {invalid_indices}")

        # For upper bounds computation, we must always use 4ti2 as Normaliz doesn't properly support
        # module generators for strict inequality problems
        # Note: self.use_4ti2 controls polymer basis computation, but upper bounds always needs 4ti2
        runner = FourTiTwoRunner()

        # We'll compute module generators for each target polymer separately
        # and combine the results (union of T_i)
        all_reactions = []

        n_on_target = len(self.on_target_indices)
        n_off_target = len(self.off_target_indices)
        on_target_list = sorted(self.on_target_indices)
        off_target_list = sorted(self.off_target_indices)

        # VARIABLE SPLITTING APPROACH TO ENFORCE CANONICAL CONSTRAINTS:
        # =============================================================
        # We transform the original n-dimensional problem into a lifted space with
        # (2 * n_on_target + n_off_target) dimensions:
        #
        # 1. On-target polymers (can be reactants or products):
        #    - Each on-target polymer p is split into two variables: p_pos and p_neg
        #    - In the original space: r[p] = p_pos - p_neg
        #    - Both p_pos >= 0 and p_neg >= 0 (enforced by Hilbert basis computation)
        #    - This allows r[p] to be positive (product), negative (reactant), or zero
        #
        # 2. Off-target polymers (can ONLY be products, not reactants):
        #    - Each off-target polymer q has a single variable q >= 0
        #    - In the original space: r[q] = q
        #    - Since q >= 0, off-target polymers can only have r[q] >= 0 (products only)
        #    - This implicitly enforces S*r >= 0 without explicitly constructing S matrix
        #
        # The lifted B matrix accounts for this variable transformation:
        B_lifted = np.zeros((self.B_matrix.shape[0], 2 * n_on_target + n_off_target), dtype=int)

        # On-target polymer columns (split into positive and negative parts)
        for i, p in enumerate(on_target_list):
            B_lifted[:, i] = self.B_matrix[:, p]  # Positive part contributes +B[p]
            B_lifted[:, n_on_target + i] = -self.B_matrix[:, p]  # Negative part contributes -B[p]

        # Off-target polymer columns (single non-negative variable)
        for i, p in enumerate(off_target_list):
            B_lifted[:, 2 * n_on_target + i] = self.B_matrix[:, p]

        # For each target polymer, compute module generators
        for target_idx in target_polymer_indices:
            # Create slice vector e_i for this target polymer
            slice_vector = np.zeros(2 * n_on_target + n_off_target, dtype=int)

            # Find position of target_idx in the lifted space
            for i, p in enumerate(off_target_list):
                if p == target_idx:
                    slice_vector[2 * n_on_target + i] = 1
                    break

            # Compute module generators for this slice
            try:
                module_gens = runner.compute_module_generators_for_slice(B_lifted, slice_vector)

                # Convert module generators from lifted space back to reaction vectors
                for h_vector in module_gens:
                    # Reconstruct the original reaction vector
                    reaction = np.zeros(n_polymers, dtype=int)

                    # On-target polymers: r[p] = h_pos[i] - h_neg[i]
                    # This can be positive (product), negative (reactant), or zero
                    for i, p in enumerate(on_target_list):
                        reaction[p] = h_vector[i] - h_vector[n_on_target + i]

                    # Off-target polymers: r[p] = h[i]
                    # Since h[i] >= 0 (from Hilbert basis), these are always non-negative
                    # This ensures off-target polymers only appear as products (canonical constraint)
                    for i, p in enumerate(off_target_list):
                        reaction[p] = h_vector[2 * n_on_target + i]

                    # Skip trivial reactions (all zeros)
                    if np.any(reaction != 0):
                        all_reactions.append(Reaction(reaction))

            except RuntimeError as e:
                print(f"Warning: Failed to compute module generators for polymer {target_idx}: {e}")
                continue

        # Remove duplicates (reactions might appear in multiple T_i)
        unique_reactions = []
        seen_vectors = set()

        for reaction in all_reactions:
            # Convert to tuple for hashing
            vector_tuple = tuple(reaction.vector)
            if vector_tuple not in seen_vectors:
                seen_vectors.add(vector_tuple)
                unique_reactions.append(reaction)

        return unique_reactions

    def check_on_target_detailed_balance(self, reactions: List[Reaction]) -> Optional[Reaction]:
        """
        Check if on-target polymers are in detailed balance.

        For reactions entirely over on-target polymers (no off-target products),
        check that they have the same number of reactants and products.

        Args:
            reactions: List of irreducible canonical reactions

        Returns:
            First violating reaction if found, None otherwise
        """
        for reaction in reactions:
            reactants, products = reaction.get_reactants_and_products()

            # Check if this reaction is entirely over on-target polymers
            all_on_target = all(idx in self.on_target_indices for idx, _ in reactants + products)

            if all_on_target and not reaction.is_balanced():
                return reaction

        return None
