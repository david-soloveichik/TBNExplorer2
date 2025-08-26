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
    """Computes irreducible canonical reactions for a TBN system."""

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

        # Split variables into positive and negative parts for on-target polymers only
        n_on_target = len(self.on_target_indices)
        n_off_target = len(self.off_target_indices)

        # Variable order: [r_on_target_pos, r_on_target_neg, r_off_target]
        # Total variables: 2 * n_on_target + n_off_target

        # Create the new B matrix for the split variables
        on_target_list = sorted(self.on_target_indices)
        off_target_list = sorted(self.off_target_indices)

        B_lifted = np.zeros((self.B_matrix.shape[0], 2 * n_on_target + n_off_target), dtype=int)

        # Fill in the lifted B matrix
        for i, p in enumerate(on_target_list):
            # Positive part
            B_lifted[:, i] = self.B_matrix[:, p]
            # Negative part (subtract)
            B_lifted[:, n_on_target + i] = -self.B_matrix[:, p]

        for i, p in enumerate(off_target_list):
            B_lifted[:, 2 * n_on_target + i] = self.B_matrix[:, p]

        # Compute Hilbert basis of { x >= 0 : B_lifted * x = 0 }
        hilbert_basis = runner.compute_hilbert_basis(B_lifted)

        if not hilbert_basis:
            return []

        # Convert back to reaction vectors
        reactions = []
        for h_vector in hilbert_basis:
            # Reconstruct the original reaction vector
            reaction = np.zeros(n_polymers, dtype=int)

            # On-target polymers: r[p] = h_pos[i] - h_neg[i]
            for i, p in enumerate(on_target_list):
                reaction[p] = h_vector[i] - h_vector[n_on_target + i]

            # Off-target polymers: r[p] = h[i]
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
        The system solved is:
        - B*r = 0 (mass conservation)
        - S*r >= 0 (canonical reactions - no off-target reactants)
        - P*r > 0 (must produce at least one target polymer)

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

        # Create P matrix: single row selecting sum of target polymers
        P_matrix = np.zeros(n_polymers, dtype=int)
        for idx in target_polymer_indices:
            P_matrix[idx] = 1

        # We need to compute Hilbert basis with strict inequality
        # System: B*r = 0, S*r >= 0, P*r > 0

        # For Normaliz, we need to transform this carefully
        # We'll use the same variable splitting approach as before but with strict inequality

        n_on_target = len(self.on_target_indices)
        n_off_target = len(self.off_target_indices)

        on_target_list = sorted(self.on_target_indices)
        off_target_list = sorted(self.off_target_indices)

        # Create lifted B matrix for split variables
        B_lifted = np.zeros((self.B_matrix.shape[0], 2 * n_on_target + n_off_target), dtype=int)

        for i, p in enumerate(on_target_list):
            B_lifted[:, i] = self.B_matrix[:, p]  # Positive part
            B_lifted[:, n_on_target + i] = -self.B_matrix[:, p]  # Negative part

        for i, p in enumerate(off_target_list):
            B_lifted[:, 2 * n_on_target + i] = self.B_matrix[:, p]

        # Create lifted S matrix (all off-target variables must be non-negative)
        # This is implicitly handled by the non-negativity constraint in Normaliz
        # But we need to ensure the transformation is correct

        # Create lifted P vector for selecting target polymers
        P_lifted = np.zeros(2 * n_on_target + n_off_target, dtype=int)
        for i, p in enumerate(off_target_list):
            if p in target_polymer_indices:
                P_lifted[2 * n_on_target + i] = 1

        # Use Normaliz - but with a filtering approach
        # The strict inequality approach isn't working, so we compute the full
        # Hilbert basis and filter for reactions that produce target polymers
        if self.use_4ti2:
            raise ValueError("Upper bound computation requires Normaliz (--use-4ti2 not supported)")

        runner = NormalizRunner()

        # Compute the regular Hilbert basis (all canonical reactions)
        hilbert_basis_raw = runner.compute_hilbert_basis(B_lifted)
        
        # Filter for reactions where P_lifted * h > 0 (i.e., at least one target polymer is produced)
        hilbert_basis = []
        for h_vector in hilbert_basis_raw:
            if np.dot(P_lifted, h_vector) > 0:
                hilbert_basis.append(h_vector)

        if not hilbert_basis:
            return []

        # Convert back to reaction vectors
        reactions = []
        for h_vector in hilbert_basis:
            # Reconstruct the original reaction vector
            reaction = np.zeros(n_polymers, dtype=int)

            # On-target polymers: r[p] = h_pos[i] - h_neg[i]
            for i, p in enumerate(on_target_list):
                reaction[p] = h_vector[i] - h_vector[n_on_target + i]

            # Off-target polymers: r[p] = h[i]
            for i, p in enumerate(off_target_list):
                reaction[p] = h_vector[2 * n_on_target + i]

            # Skip trivial reactions (all zeros) - shouldn't happen with strict inequality
            if np.any(reaction != 0):
                reactions.append(Reaction(reaction))

        return reactions

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
