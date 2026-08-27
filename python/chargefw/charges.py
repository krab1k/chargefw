"""Owned charge assignment values."""

from __future__ import annotations

from collections.abc import Hashable
from dataclasses import dataclass

import numpy as np

from .core import SourceIdentity


def _readonly_values(values: np.ndarray) -> np.ndarray:
    """Copy values into an array backed by immutable Python bytes."""
    array = np.asarray(values, dtype=np.float64, order="C")
    return np.frombuffer(array.tobytes(), dtype=np.float64).reshape(array.shape)


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
        object.__setattr__(self, "values", _readonly_values(self.values))

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
        result = object.__new__(cls)
        object.__setattr__(result, "values", _readonly_values(values))
        object.__setattr__(result, "molecule_index", molecule_index)
        object.__setattr__(result, "conformer_index", conformer_index)
        object.__setattr__(result, "source", source)
        object.__setattr__(result, "atom_ids", atom_ids)
        object.__setattr__(result, "conformer_id", conformer_id)
        return result
