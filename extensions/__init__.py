"""Extensions package for TBNExplorer2.

This package contains additional functionality that builds on the core TBN framework,
including:
- Enumeration of irreducible canonical reactions
- IBOT (Iterative Balancing of Off-Target) algorithm
"""

from .canonical_reactions import CanonicalReactionsComputer
from .ibot import IBOTAlgorithm

__all__ = ["CanonicalReactionsComputer", "IBOTAlgorithm"]
