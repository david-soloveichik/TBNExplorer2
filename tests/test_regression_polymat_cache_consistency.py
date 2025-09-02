from pathlib import Path

import numpy as np

from tbnexplorer2.model import BindingSite, Monomer, TBN
from tbnexplorer2.polymer_basis import Polymer, PolymerBasisComputer
from tbnexplorer2.polymat_io import PolymatReader


class FakeNupackRunner:
    """Deterministic concentration runner emulating NUPACK behavior for tests.

    Concentrations are proportional to exp(-G), where G is the polymer free energy
    computed with the provided deltaG and temperature.
    """

    def __init__(self, temperature: float = 37.0):
        self.temperature = temperature

    def check_nupack_available(self) -> bool:
        return True

    def compute_equilibrium_concentrations(self, polymers, tbn, output_dir=None, deltaG=None, temperature=37.0):
        energies = np.array([p.compute_free_energy(deltaG, temperature) for p in polymers], dtype=float)
        # Softmax-like but without normalization (relative values only matter for sorting)
        return np.exp(-energies)


def test_concentrations_identical_with_and_without_cached_polymat(tmp_path):
    """Regression: concentrations should not depend on presence of .tbnpolymat.

    This mimics running with --use-nupack-concentrations and --deltaG by supplying
    a fake NUPACK-like runner, and verifies results are identical across runs
    whether the basis is freshly computed or loaded from an existing .tbnpolymat.
    """
    # Build simple TBN with units and concentrations
    # Monomer1: a, Monomer2: a*
    m1 = Monomer(name="M1", binding_sites=[BindingSite("a", False)], concentration=100.0, original_line="M1: a, 100")
    m2 = Monomer(name="M2", binding_sites=[BindingSite("a", True)], concentration=100.0, original_line="M2: a*, 100")
    tbn = TBN([m1, m2], {"a": 0}, concentration_units="nM")

    # Define a few polymers (vector counts)
    polymers = [
        Polymer(np.array([1, 1]), [m1, m2], tbn),  # dimer
        Polymer(np.array([2, 0]), [m1, m2], tbn),  # two M1
        Polymer(np.array([0, 2]), [m1, m2], tbn),  # two M2
    ]

    computer = PolymerBasisComputer(tbn)
    out_file = Path(tmp_path) / "example.tbnpolymat"

    # Runner that uses current deltaG to compute energies deterministically
    runner = FakeNupackRunner(temperature=37.0)
    deltaG = [-2.0, 5.0, 3.0]

    # First save (no cache yet)
    computer.save_tbnpolymat(
        polymers,
        str(out_file),
        compute_free_energies=True,
        compute_concentrations=True,
        concentration_runner=runner,
        verbose=False,
        deltaG=deltaG,
        temperature=37.0,
    )

    r1 = PolymatReader(str(out_file)).read()
    c1 = r1.concentrations.copy()

    # Now load from cache and save again to the same file
    polymers_cached = computer.load_cached_polymer_basis(str(out_file))
    assert polymers_cached is not None

    computer.save_tbnpolymat(
        polymers_cached,
        str(out_file),
        compute_free_energies=True,
        compute_concentrations=True,
        concentration_runner=runner,
        verbose=False,
        deltaG=deltaG,
        temperature=37.0,
    )

    r2 = PolymatReader(str(out_file)).read()
    c2 = r2.concentrations.copy()

    # Concentrations should be identical across runs
    assert np.allclose(c1, c2)
