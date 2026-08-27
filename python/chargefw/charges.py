"""Owned charge assignment values."""

from __future__ import annotations

from collections.abc import Hashable
from dataclasses import dataclass

import numpy as np

from .core import SourceIdentity


@dataclass(frozen=True, slots=True)
class ChargeAssignment:
    """Source-mapped, owned charges for one molecule or conformer."""

    values: np.ndarray
    molecule_index: int
    conformer_index: int | None
    source: SourceIdentity
    atom_ids: tuple[Hashable, ...]
    conformer_id: Hashable | None

    def __post_init__(self) -> None:
        values = np.array(self.values, dtype=np.float64, order="C", copy=True)
        values.setflags(write=False)
        object.__setattr__(self, "values", values)

    @classmethod
    def _from_native_values(
        cls,
        values: np.ndarray,
        *,
        molecule_index: int,
        conformer_index: int | None,
        source: SourceIdentity,
        atom_ids: tuple[Hashable, ...],
        conformer_id: Hashable | None,
    ) -> "ChargeAssignment":
        """Build an assignment from an extension-owned, one-dimensional NumPy array."""
        if values.dtype != np.float64 or values.ndim != 1 or not values.flags.c_contiguous:
            raise RuntimeError("native charge values must be a C-contiguous float64 vector")
        values.setflags(write=False)
        result = object.__new__(cls)
        object.__setattr__(result, "values", values)
        object.__setattr__(result, "molecule_index", molecule_index)
        object.__setattr__(result, "conformer_index", conformer_index)
        object.__setattr__(result, "source", source)
        object.__setattr__(result, "atom_ids", atom_ids)
        object.__setattr__(result, "conformer_id", conformer_id)
        return result
