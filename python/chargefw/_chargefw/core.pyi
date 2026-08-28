from collections.abc import Sequence

import numpy as np
from numpy.typing import NDArray

class _NativeMolecule:
    @property
    def name(self) -> str: ...

    @property
    def atom_count(self) -> int: ...

    @property
    def bond_count(self) -> int: ...

    @property
    def conformer_count(self) -> int: ...


def _make_molecule(
    atomic_numbers: NDArray[np.int64],
    formal_charges: NDArray[np.int64],
    bonds: NDArray[np.int64],
    coordinates: NDArray[np.float64],
    atom_names: Sequence[str],
    conformer_names: Sequence[str],
    name: str,
) -> _NativeMolecule: ...
