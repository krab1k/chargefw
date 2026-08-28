"""Owned charge assignment values."""

from __future__ import annotations

from collections.abc import Hashable, Iterable
from dataclasses import dataclass
from operator import index as as_index
from typing import Any

import numpy as np

from ._arrays import immutable_array
from .core import SourceIdentity


def _charge_values(values: Any) -> np.ndarray:
    try:
        array = np.asarray(values)
    except Exception as error:
        raise TypeError("values must be a real numeric vector") from error
    if array.ndim != 1:
        raise ValueError("values must be a one-dimensional vector")
    if array.dtype.kind not in "iuf" or array.dtype.kind == "b":
        raise TypeError("values must contain real numbers")
    try:
        normalized = np.asarray(array, dtype=np.float64, order="C")
    except (OverflowError, TypeError, ValueError) as error:
        raise ValueError("values contain values outside the supported range") from error
    if not np.all(np.isfinite(normalized)):
        raise ValueError("values must contain only finite charges")
    return immutable_array(normalized)


def _nonnegative_index(value: Any, name: str) -> int:
    if isinstance(value, (bool, np.bool_)):
        raise TypeError(f"{name} must be an integer")
    try:
        result = as_index(value)
    except TypeError as error:
        raise TypeError(f"{name} must be an integer") from error
    if result < 0:
        raise ValueError(f"{name} must be non-negative")
    return result


def _hashable_ids(values: Iterable[Hashable], count: int) -> tuple[Hashable, ...]:
    if isinstance(values, (str, bytes)):
        raise TypeError("atom_ids must be an iterable of hashable values")
    try:
        result = tuple(values)
    except TypeError as error:
        raise TypeError("atom_ids must be an iterable of hashable values") from error
    if len(result) != count:
        raise ValueError("atom_ids must have the same length as values")
    for value in result:
        try:
            hash(value)
        except TypeError as error:
            raise ValueError("atom_ids must contain only hashable values") from error
    return result


def _validate_hashable(value: Hashable | None, name: str) -> None:
    if value is None:
        return
    try:
        hash(value)
    except TypeError as error:
        raise ValueError(f"{name} must be hashable") from error


@dataclass(frozen=True, slots=True, eq=False)
class ChargeAssignment:
    """Source-mapped, owned charges for one molecule or conformer."""

    values: np.ndarray
    molecule_index: int
    conformer_index: int | None
    source: SourceIdentity
    atom_ids: tuple[Hashable, ...]
    conformer_id: Hashable | None

    def __post_init__(self) -> None:
        values = _charge_values(self.values)
        molecule_index = _nonnegative_index(self.molecule_index, "molecule_index")
        conformer_index = (
            None
            if self.conformer_index is None
            else _nonnegative_index(self.conformer_index, "conformer_index")
        )
        if not isinstance(self.source, SourceIdentity):
            raise TypeError("source must be a SourceIdentity")
        atom_ids = _hashable_ids(self.atom_ids, len(values))
        _validate_hashable(self.conformer_id, "conformer_id")
        if conformer_index is None and self.conformer_id is not None:
            raise ValueError("conformer_id requires a conformer_index")

        object.__setattr__(self, "values", values)
        object.__setattr__(self, "molecule_index", molecule_index)
        object.__setattr__(self, "conformer_index", conformer_index)
        object.__setattr__(self, "atom_ids", atom_ids)

    def __eq__(self, other: object) -> bool:
        if not isinstance(other, ChargeAssignment):
            return NotImplemented
        return (
            np.array_equal(self.values, other.values)
            and self.molecule_index == other.molecule_index
            and self.conformer_index == other.conformer_index
            and self.source == other.source
            and self.atom_ids == other.atom_ids
            and self.conformer_id == other.conformer_id
        )

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
        if len(values) != len(atom_ids):
            raise RuntimeError("native charge count does not match the source atom mapping")
        if not np.all(np.isfinite(values)):
            raise RuntimeError("native charge values must be finite")
        result = object.__new__(cls)
        object.__setattr__(result, "values", immutable_array(values))
        object.__setattr__(result, "molecule_index", molecule_index)
        object.__setattr__(result, "conformer_index", conformer_index)
        object.__setattr__(result, "source", source)
        object.__setattr__(result, "atom_ids", atom_ids)
        object.__setattr__(result, "conformer_id", conformer_id)
        return result
