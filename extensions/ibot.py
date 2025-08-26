"""IBOT (Iterative Balancing of Off-Target) algorithm implementation.

This module implements the IBOT algorithm for assigning concentration exponents
to off-target polymers such that they are in detailed balance.
"""

from dataclasses import dataclass
from pathlib import Path
from typing import Dict, List, Set

import numpy as np

from tbnexplorer2.model import TBN
from tbnexplorer2.tbnpolys_io import TbnpolysWriter
from tbnexplorer2.units import to_molar

from .canonical_reactions import Reaction


@dataclass
class ReactionMetrics:
    """Metrics for a reaction in the IBOT algorithm."""

    novelty: int  # Number of unassigned off-target polymers
    imbalance: float  # Concentration exponent imbalance
    ratio: float  # imbalance/novelty ratio


@dataclass
class IterationInfo:
    """Information about a single IBOT iteration."""

    iteration: int
    mu_min: float
    reactions: List[Reaction]  # Reactions that were part of R in this iteration
    assigned_polymers: Set[int]  # Polymer indices assigned in this iteration


class IBOTAlgorithm:
    """Implementation of the Iterative Balancing of Off-Target algorithm."""

    def __init__(self, tbn: TBN, polymers: List[np.ndarray], on_target_indices: Set[int], reactions: List[Reaction]):
        """
        Initialize IBOT algorithm.

        Args:
            tbn: The TBN model
            polymers: List of all polymers in the basis
            on_target_indices: Set of indices of on-target polymers
            reactions: List of irreducible canonical reactions
        """
        self.tbn = tbn
        self.polymers = polymers
        self.on_target_indices = on_target_indices
        self.off_target_indices = set(range(len(polymers))) - on_target_indices
        self.reactions = reactions

        # Initialize concentration exponents
        self.mu = np.zeros(len(polymers))
        # On-target polymers have μ = 1
        for idx in on_target_indices:
            self.mu[idx] = 1.0

        # Track which off-target polymers have been assigned
        self.unassigned_off_target = self.off_target_indices.copy()

        # Track iteration information for reactions output
        self.iteration_info = []

    def compute_reaction_metrics(self, reaction: Reaction) -> ReactionMetrics:
        """
        Compute novelty and imbalance for a reaction.

        Args:
            reaction: The reaction to analyze

        Returns:
            ReactionMetrics object with novelty, imbalance, and ratio
        """
        novelty = 0
        imbalance = 0.0

        for i, count in enumerate(reaction.vector):
            if count != 0:
                # Check if this is an unassigned off-target polymer
                if i in self.unassigned_off_target:
                    novelty += 1

                # Compute contribution to imbalance
                if count < 0:  # Reactant
                    imbalance += abs(count) * self.mu[i]
                else:  # Product
                    imbalance -= count * self.mu[i]

        # Compute ratio (handle division by zero)
        ratio = imbalance / novelty if novelty > 0 else float("inf")

        return ReactionMetrics(novelty, imbalance, ratio)

    def run(self) -> Dict[int, float]:
        """
        Run the IBOT algorithm to assign concentration exponents.

        Returns:
            Dictionary mapping polymer index to concentration exponent
        """
        iteration = 0

        while self.unassigned_off_target:
            iteration += 1

            # Compute metrics for all reactions with novelty > 0
            active_reactions = []
            for reaction in self.reactions:
                metrics = self.compute_reaction_metrics(reaction)
                if metrics.novelty > 0:
                    active_reactions.append((reaction, metrics))

            if not active_reactions:
                # No more reactions with unassigned polymers
                break

            # Find minimum imbalance-novelty ratio
            min_ratio = min(metrics.ratio for _, metrics in active_reactions)

            # Find all reactions with this minimum ratio
            min_reactions = [
                reaction for reaction, metrics in active_reactions if abs(metrics.ratio - min_ratio) < 1e-10
            ]

            # Collect all unassigned off-target polymers appearing in these reactions
            polymers_to_assign = set()
            for reaction in min_reactions:
                for i, count in enumerate(reaction.vector):
                    if count != 0 and i in self.unassigned_off_target:
                        polymers_to_assign.add(i)

            # Store iteration information
            iter_info = IterationInfo(
                iteration=iteration,
                mu_min=min_ratio,
                reactions=min_reactions.copy(),
                assigned_polymers=polymers_to_assign.copy(),
            )
            self.iteration_info.append(iter_info)

            # Assign concentration exponent to these polymers
            for p in polymers_to_assign:
                self.mu[p] = min_ratio
                self.unassigned_off_target.remove(p)

            print(f"IBOT iteration {iteration}: Assigned μ={min_ratio:.6f} to {len(polymers_to_assign)} polymers")

        # Remove polymers that were never assigned (μ = 0)
        # These are off-target polymers that cannot be generated from on-target ones
        assigned_polymers = {}
        for i in range(len(self.polymers)):
            # Keep all on-target polymers and off-target polymers with μ > 0
            if i in self.on_target_indices or self.mu[i] > 0:
                assigned_polymers[i] = self.mu[i]

        return assigned_polymers

    def generate_tbnpolys_output(self, output_file: Path):
        """
        Generate .tbnpolys file with concentration exponents.

        Args:
            output_file: Path to output .tbnpolys file
        """
        writer = TbnpolysWriter(self.tbn)

        # Separate on-target and off-target polymers
        on_target_polymers = []
        on_target_mus = []
        off_target_polymers = []
        off_target_mus = []

        for i, polymer in enumerate(self.polymers):
            # Skip unassigned off-target polymers (μ = 0)
            if i not in self.on_target_indices and self.mu[i] == 0:
                continue

            if i in self.on_target_indices:
                on_target_polymers.append(polymer)
                on_target_mus.append(self.mu[i])
            else:
                off_target_polymers.append(polymer)
                off_target_mus.append(self.mu[i])

        # Sort off-target polymers by concentration exponent (ascending)
        if off_target_polymers:
            sorted_indices = np.argsort(off_target_mus)
            off_target_polymers = [off_target_polymers[i] for i in sorted_indices]
            off_target_mus = [off_target_mus[i] for i in sorted_indices]

        # Write file content
        lines = []

        # Header
        lines.append("# IBOT Results - Concentration Exponents")
        lines.append(f"# Total polymers: {len(self.polymers)}")
        lines.append(f"# On-target polymers: {len(on_target_polymers)}")
        lines.append(f"# Off-target polymers: {len(off_target_polymers)}")
        lines.append("")

        # On-target polymers
        lines.append("# === ON-TARGET POLYMERS ===")
        lines.append("")

        for polymer, mu_val in zip(on_target_polymers, on_target_mus):
            polymer_lines = self._format_polymer_with_mu(polymer, mu_val, writer)
            lines.extend(polymer_lines)
            lines.append("")

        # Off-target polymers
        if off_target_polymers:
            lines.append("# === OFF-TARGET POLYMERS ===")
            lines.append("# (sorted by concentration exponent)")
            lines.append("")

            for polymer, mu_val in zip(off_target_polymers, off_target_mus):
                polymer_lines = self._format_polymer_with_mu(polymer, mu_val, writer)
                lines.extend(polymer_lines)
                lines.append("")

        # Remove trailing empty line
        if lines and lines[-1] == "":
            lines.pop()

        # Write to file
        with open(output_file, "w") as f:
            f.write("\n".join(lines))

        print(f"Saved IBOT results to {output_file}")

    def generate_reactions_output(self, output_file: Path):
        """
        Generate text file showing irreducible canonical reactions.

        Reactions are ordered by the iteration of the IBOT algorithm in which they
        were part of R. The μ_min value and polymers assigned in each iteration are indicated.

        Args:
            output_file: Path to output text file
        """
        lines = []

        # Header
        lines.append("# Irreducible Canonical Reactions from IBOT Algorithm")
        lines.append(f"# Total reactions: {len(self.reactions)}")
        lines.append(f"# Total iterations: {len(self.iteration_info)}")
        lines.append("#")
        lines.append("# Notation:")
        lines.append("#   - Polymers are shown in brackets: {monomer1; monomer2; ...}")
        lines.append("#   - Monomer multiplicities shown as prefix: {2 monomer1; monomer2}")
        lines.append("#   - Polymers marked with ^ were assigned μ in that iteration")
        lines.append("=" * 80)
        lines.append("")

        # Process each iteration
        for iter_info in self.iteration_info:
            lines.append(f"## Iteration {iter_info.iteration}")
            lines.append(f"## μ_min = {iter_info.mu_min:.6f}")
            lines.append(f"## Number of reactions in R: {len(iter_info.reactions)}")
            lines.append(f"## Polymers assigned μ in this iteration: {len(iter_info.assigned_polymers)}")
            lines.append("")

            # Format each reaction in this iteration
            for reaction in iter_info.reactions:
                reaction_str = self._format_reaction_with_assignments(reaction, iter_info.assigned_polymers)
                lines.append(reaction_str)

            lines.append("")
            lines.append("-" * 40)
            lines.append("")

        # Remove trailing separators
        while lines and lines[-1] in ["", "-" * 40]:
            lines.pop()

        # Write to file
        with open(output_file, "w") as f:
            f.write("\n".join(lines))

        print(f"Saved canonical reactions output to {output_file}")

    def _format_reaction_with_assignments(self, reaction: Reaction, assigned_polymers: Set[int]) -> str:
        """
        Format a reaction with indicators for newly assigned polymers.

        Args:
            reaction: The reaction to format
            assigned_polymers: Set of polymer indices assigned in this iteration

        Returns:
            Formatted reaction string
        """
        reactants, products = reaction.get_reactants_and_products()

        def format_side(polymers):
            terms = []
            for idx, mult in polymers:
                # Create polymer representation
                polymer_str = self._get_polymer_representation(idx)

                # Mark if this polymer was assigned in this iteration
                if idx in assigned_polymers:
                    polymer_str += "^"  # Mark with caret to indicate newly assigned μ

                if mult == 1:
                    terms.append(polymer_str)
                else:
                    terms.append(f"{mult} {polymer_str}")
            return " + ".join(terms) if terms else "0"

        return f"{format_side(reactants)} -> {format_side(products)}"

    def _get_polymer_representation(self, polymer_idx: int) -> str:
        """
        Get a string representation of a polymer in .tbnpolys notation.

        Args:
            polymer_idx: Index of the polymer

        Returns:
            String representation in format {monomer1; monomer2; ...}
        """
        polymer = self.polymers[polymer_idx]

        # Collect all monomers in the polymer
        monomer_specs = []
        for monomer_idx, count in enumerate(polymer):
            if count > 0:
                monomer = self.tbn.monomers[monomer_idx]

                # Get monomer representation (name or binding sites)
                monomer_spec = monomer.name or monomer.get_binding_sites_str()

                # Add multiplicity prefix if count > 1
                if count > 1:
                    monomer_spec = f"{count} {monomer_spec}"

                monomer_specs.append(monomer_spec)

        # Join with semicolons and wrap in brackets
        return "{" + "; ".join(monomer_specs) + "}"

    def _format_polymer_with_mu(self, polymer: np.ndarray, mu_val: float, writer: TbnpolysWriter) -> List[str]:
        """
        Format a polymer with its concentration exponent.

        Args:
            polymer: Polymer monomer counts
            mu_val: Concentration exponent value
            writer: TbnpolysWriter instance for formatting

        Returns:
            List of lines representing the polymer
        """
        lines = writer._format_single_polymer(polymer)
        lines.append(f"# μ: {mu_val:.6f}")
        return lines

    def generate_tbn_output(self, output_file: Path, c: float, units: str):
        """
        Generate .tbn file with computed monomer concentrations.

        Each monomer i is assigned concentration sum over all polymers p of:
        p[i] * ((c'/rho_H2O)^μ(p)) * rho_H2O where p[i] is the count of monomer i in polymer p,
        c' is c converted to Molar, and rho_H2O = 55.14 M is the density of water at 37C.

        Args:
            output_file: Path to output .tbn file
            c: Base concentration value (in specified units)
            units: Concentration units (e.g., 'nM', 'uM')
        """
        # Convert c to Molar for computation
        c_molar = to_molar(c, units)

        # Water density at 37C in Molar
        rho_h2o = 55.14

        # Compute monomer concentrations in Molar
        monomer_concentrations = np.zeros(len(self.tbn.monomers))

        for p_idx, polymer in enumerate(self.polymers):
            # Skip unassigned off-target polymers
            if p_idx not in self.on_target_indices and self.mu[p_idx] == 0:
                continue

            mu_p = self.mu[p_idx]
            # Use mole fraction logic: ((c'/rho_H2O)^μ(p)) * rho_H2O
            concentration_factor = ((c_molar / rho_h2o) ** mu_p) * rho_h2o

            for m_idx, count in enumerate(polymer):
                if count > 0:
                    monomer_concentrations[m_idx] += count * concentration_factor

        # Convert concentrations back to specified units
        from tbnexplorer2.units import from_molar

        monomer_concentrations = from_molar(monomer_concentrations, units)

        # Generate .tbn file content
        lines = []

        # Add UNITS line
        lines.append(f"\\UNITS: {units}")
        lines.append("")

        # Add monomers with concentrations
        for m_idx, monomer in enumerate(self.tbn.monomers):
            concentration = monomer_concentrations[m_idx]

            # Format monomer line
            if monomer.name:
                # Use the original format with name
                if ":" in monomer.original_line:
                    # Name before binding sites
                    monomer_line = f"{monomer.name}: {monomer.get_binding_sites_str()}"
                else:
                    # Name after binding sites
                    monomer_line = f"{monomer.get_binding_sites_str()} >{monomer.name}"
            else:
                # No name, just binding sites
                monomer_line = monomer.get_binding_sites_str()

            # Add concentration
            monomer_line += f", {concentration:.6g}"
            lines.append(monomer_line)

        # Write to file
        with open(output_file, "w") as f:
            f.write("\n".join(lines))

        print(f"Generated .tbn file with concentrations at {output_file}")
        print(f"Base concentration c = {c} {units}")
