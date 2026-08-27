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
    source_atom_ids: tuple[Hashable, ...]
    source_conformer_id: Hashable | None

    def __post_init__(self) -> None:
        values = np.array(self.values, dtype=np.float64, order="C", copy=True)
        values.setflags(write=False)
        object.__setattr__(self, "values", values)

    @property
    def atom_ids(self) -> tuple[Hashable, ...]:
        return self.source_atom_ids
