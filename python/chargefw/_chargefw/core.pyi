from collections.abc import Sequence

import numpy as np
from numpy.typing import NDArray

class _NativeMolecule: ...

def _make_molecule(
    atomic_numbers: NDArray[np.int64],
    formal_charges: NDArray[np.int64],
    bonds: NDArray[np.int64],
    coordinates: NDArray[np.float64],
    atom_names: Sequence[str],
    conformer_names: Sequence[str],
    name: str,
) -> _NativeMolecule: ...
