"""ChargeFW Python API."""

from __future__ import annotations

from collections.abc import Hashable, Iterable, Sequence
from dataclasses import dataclass
import numbers
from operator import index as as_index
from typing import Any

import numpy as np

from ._chargefw import version as _native_version
from . import _chargefw

__version__ = _native_version()


@dataclass(frozen=True, slots=True)
class SourceIdentity:
    """Identity of the source record represented by a molecule."""

    source: str = ""
    record_index: int = 0
    record_id: Hashable | None = None

    def __post_init__(self) -> None:
        if not isinstance(self.source, str):
            raise TypeError("source must be a string")
        try:
            normalized_record_index = as_index(self.record_index)
        except TypeError as error:
            raise TypeError("record_index must be an integer") from error
        if normalized_record_index < 0:
            raise ValueError("record_index must be non-negative")
        if self.record_id is not None:
            try:
                hash(self.record_id)
            except TypeError as error:
                raise ValueError("record_id must be hashable") from error
        object.__setattr__(self, "record_index", normalized_record_index)

    @property
    def source_name(self) -> str:
        return self.source


def _as_integer_array(value: Any, field: str, ndim: int) -> np.ndarray:
    try:
        array = np.asarray(value)
    except Exception as error:
        raise TypeError(f"{field} must be an integer array") from error

    if array.ndim != ndim:
        raise ValueError(f"{field} must have rank {ndim}, got {array.ndim}")
    if array.dtype.kind == "b":
        raise TypeError(f"{field} must contain integers, not booleans")
    if array.size == 0 and array.dtype.kind == "f":
        return np.empty(array.shape, dtype=np.int64)
    if array.dtype.kind not in "iu":
        if array.dtype.kind == "O" and all(
            isinstance(item, numbers.Integral) and not isinstance(item, (bool, np.bool_))
            for item in array.flat
        ):
            raise ValueError(f"{field} contains values outside the supported integer range")
        raise TypeError(f"{field} must contain integers")

    if array.dtype.kind == "u" and np.any(array > np.iinfo(np.int64).max):
        raise ValueError(f"{field} contains values outside the supported integer range")
    try:
        return np.array(array, dtype=np.int64, order="C", copy=True)
    except (OverflowError, ValueError) as error:
        raise ValueError(f"{field} contains values outside the supported integer range") from error


def _as_coordinates(value: Any, atom_count: int) -> tuple[np.ndarray, bool]:
    if value is None:
        return np.empty((0, atom_count, 3), dtype=np.float64), False
    try:
        array = np.asarray(value)
    except Exception as error:
        raise TypeError("coordinates must be a numeric array") from error
    if array.ndim == 2:
        if array.shape != (atom_count, 3):
            raise ValueError(f"coordinates must have shape ({atom_count}, 3)")
        one_conformer = True
        canonical_shape = (1, atom_count, 3)
    elif array.ndim == 3:
        if array.shape[1:] != (atom_count, 3):
            raise ValueError(f"coordinates must have shape (C, {atom_count}, 3)")
        one_conformer = False
        canonical_shape = array.shape
    else:
        raise ValueError("coordinates must have rank 2 or 3")
    if array.dtype.kind not in "iuf" or array.dtype.kind == "b":
        raise TypeError("coordinates must contain real numbers")
    try:
        canonical = np.array(array, dtype=np.float64, order="C", copy=True).reshape(canonical_shape)
    except (OverflowError, ValueError) as error:
        raise ValueError("coordinates contain values outside the supported range") from error
    if not np.all(np.isfinite(canonical)):
        raise ValueError("coordinates must contain only finite values")
    public = canonical.reshape((atom_count, 3)) if one_conformer else canonical
    return public, one_conformer


def _as_names(value: Sequence[str] | None, count: int, field: str) -> tuple[str, ...]:
    if value is None:
        return ("",) * count
    if isinstance(value, (str, bytes)):
        raise TypeError(f"{field} must be a sequence of strings")
    try:
        result = tuple(value)
    except TypeError as error:
        raise TypeError(f"{field} must be a sequence of strings") from error
    if len(result) != count:
        raise ValueError(f"{field} must contain {count} values")
    if not all(isinstance(item, str) for item in result):
        raise TypeError(f"{field} must contain only strings")
    return result


def _as_ids(value: Iterable[Hashable] | None, count: int, field: str) -> tuple[Hashable, ...]:
    if value is None:
        return tuple(range(count))
    if isinstance(value, (str, bytes)):
        raise TypeError(f"{field} must be a sequence of hashable values")
    try:
        result = tuple(value)
    except TypeError as error:
        raise TypeError(f"{field} must be a sequence of hashable values") from error
    if len(result) != count:
        raise ValueError(f"{field} must contain {count} values")
    for item in result:
        try:
            hash(item)
        except TypeError as error:
            raise ValueError(f"{field} must contain only hashable values") from error
    return result


def _readonly_copy(array: np.ndarray) -> np.ndarray:
    result = array.copy()
    result.setflags(write=False)
    return result


class Molecule:
    """Immutable, owned toolkit-neutral molecular data."""

    __slots__ = (
        "_atomic_numbers",
        "_formal_charges",
        "_bonds",
        "_coordinates",
        "_native_coordinates",
        "_coordinates_were_2d",
        "_name",
        "_atom_names",
        "_conformer_names",
        "_source",
        "_atom_ids",
        "_conformer_ids",
        "_native",
    )

    def __init__(
        self,
        atomic_numbers: Any,
        formal_charges: Any = None,
        bonds: Any = None,
        coordinates: Any = None,
        name: str | None = None,
        atom_names: Sequence[str] | None = None,
        conformer_names: Sequence[str] | None = None,
        source: SourceIdentity | None = None,
        source_name: str | None = None,
        record_index: int = 0,
        record_id: Hashable | None = None,
        atom_ids: Iterable[Hashable] | None = None,
        conformer_ids: Iterable[Hashable] | None = None,
        source_atom_ids: Iterable[Hashable] | None = None,
        source_conformer_ids: Iterable[Hashable] | None = None,
    ) -> None:
        if name is not None and not isinstance(name, str):
            raise TypeError("name must be a string or None")
        atomic = _as_integer_array(atomic_numbers, "atomic_numbers", 1)
        if np.any((atomic < 1) | (atomic > 118)):
            raise ValueError("atomic_numbers must contain values in range 1..118")
        atom_count = len(atomic)

        if formal_charges is None:
            formal = np.zeros(atom_count, dtype=np.int64)
        else:
            formal = _as_integer_array(formal_charges, "formal_charges", 1)
            if len(formal) != atom_count:
                raise ValueError("formal_charges must have the same length as atomic_numbers")
        int_limits = np.iinfo(np.int32)
        if np.any((formal < int_limits.min) | (formal > int_limits.max)):
            raise ValueError("formal_charges contain values outside the native integer range")

        if bonds is None:
            bond_array = np.empty((0, 3), dtype=np.int64)
        else:
            bond_array = _as_integer_array(bonds, "bonds", 2)
            if bond_array.shape[1] != 3:
                raise ValueError("bonds must have shape (B, 3)")
        seen_bonds: set[tuple[int, int]] = set()
        for first, second, order in bond_array:
            if first < 0 or second < 0 or first >= atom_count or second >= atom_count:
                raise ValueError("bond atom indices must be within the molecule")
            if first == second:
                raise ValueError("bond endpoints must refer to different atoms")
            if order not in (1, 2, 3):
                raise ValueError("bond orders must be 1, 2, or 3")
            key = (min(int(first), int(second)), max(int(first), int(second)))
            if key in seen_bonds:
                raise ValueError("molecule contains duplicate bonds")
            seen_bonds.add(key)

        public_coordinates, coordinates_were_2d = _as_coordinates(coordinates, atom_count)
        native_coordinates = public_coordinates
        if coordinates_were_2d:
            native_coordinates = public_coordinates.reshape((1, atom_count, 3))
        atom_name_values = _as_names(atom_names, atom_count, "atom_names")
        conformer_count = native_coordinates.shape[0]
        conformer_name_values = _as_names(conformer_names, conformer_count, "conformer_names")

        if source is not None:
            if not isinstance(source, SourceIdentity):
                raise TypeError("source must be a SourceIdentity")
            if source_name is not None or record_index != 0 or record_id is not None:
                raise ValueError("source cannot be combined with source identity fields")
            source_identity = source
        else:
            if source_name is not None and not isinstance(source_name, str):
                raise TypeError("source_name must be a string or None")
            try:
                normalized_record_index = as_index(record_index)
            except TypeError as error:
                raise TypeError("record_index must be an integer") from error
            if normalized_record_index < 0:
                raise ValueError("record_index must be non-negative")
            if record_id is not None:
                try:
                    hash(record_id)
                except TypeError as error:
                    raise ValueError("record_id must be hashable") from error
            source_identity = SourceIdentity(source_name or "", normalized_record_index, record_id)

        if atom_ids is not None and source_atom_ids is not None:
            raise ValueError("atom_ids and source_atom_ids are aliases and cannot both be supplied")
        if conformer_ids is not None and source_conformer_ids is not None:
            raise ValueError(
                "conformer_ids and source_conformer_ids are aliases and cannot both be supplied"
            )
        atom_id_values = _as_ids(
            atom_ids if atom_ids is not None else source_atom_ids, atom_count, "atom_ids"
        )
        conformer_id_values = _as_ids(
            conformer_ids if conformer_ids is not None else source_conformer_ids,
            conformer_count,
            "conformer_ids",
        )

        canonical_coordinates = native_coordinates
        object.__setattr__(self, "_atomic_numbers", atomic)
        object.__setattr__(self, "_formal_charges", formal)
        object.__setattr__(self, "_bonds", bond_array)
        object.__setattr__(self, "_coordinates", public_coordinates)
        object.__setattr__(self, "_native_coordinates", canonical_coordinates)
        object.__setattr__(self, "_coordinates_were_2d", coordinates_were_2d)
        object.__setattr__(self, "_name", name or "")
        object.__setattr__(self, "_atom_names", atom_name_values)
        object.__setattr__(self, "_conformer_names", conformer_name_values)
        object.__setattr__(self, "_source", source_identity)
        object.__setattr__(self, "_atom_ids", atom_id_values)
        object.__setattr__(self, "_conformer_ids", conformer_id_values)
        object.__setattr__(
            self,
            "_native",
            _chargefw._make_molecule(
                atomic,
                formal,
                bond_array,
                canonical_coordinates,
                atom_name_values,
                conformer_name_values,
                name or "",
            ),
        )

    def __setattr__(self, name: str, value: Any) -> None:
        raise AttributeError(f"{type(self).__name__} is immutable")

    @property
    def atomic_numbers(self) -> np.ndarray:
        return _readonly_copy(self._atomic_numbers)

    @property
    def formal_charges(self) -> np.ndarray:
        return _readonly_copy(self._formal_charges)

    @property
    def bonds(self) -> np.ndarray:
        return _readonly_copy(self._bonds)

    @property
    def coordinates(self) -> np.ndarray:
        return _readonly_copy(self._coordinates)

    @property
    def name(self) -> str:
        return self._name

    @property
    def atom_names(self) -> tuple[str, ...]:
        return self._atom_names

    @property
    def conformer_names(self) -> tuple[str, ...]:
        return self._conformer_names

    @property
    def source(self) -> SourceIdentity:
        return self._source

    @property
    def source_name(self) -> str:
        return self._source.source

    @property
    def record_index(self) -> int:
        return self._source.record_index

    @property
    def record_id(self) -> Hashable | None:
        return self._source.record_id

    @property
    def atom_ids(self) -> tuple[Hashable, ...]:
        return self._atom_ids

    @property
    def conformer_ids(self) -> tuple[Hashable, ...]:
        return self._conformer_ids

    @property
    def source_atom_ids(self) -> tuple[Hashable, ...]:
        return self._atom_ids

    @property
    def source_conformer_ids(self) -> tuple[Hashable, ...]:
        return self._conformer_ids

    @property
    def atom_count(self) -> int:
        return len(self._atomic_numbers)

    @property
    def bond_count(self) -> int:
        return len(self._bonds)

    @property
    def conformer_count(self) -> int:
        return self._native_coordinates.shape[0]

    @property
    def has_coordinates(self) -> bool:
        return self.conformer_count != 0


class MoleculeCollection(Sequence[Molecule]):
    """Immutable source-ordered collection of molecules."""

    __slots__ = ("_molecules", "_name", "_native")

    def __init__(self, molecules: Iterable[Molecule], name: str | None = None) -> None:
        if isinstance(molecules, (str, bytes, Molecule)):
            raise TypeError("molecules must be an iterable of Molecule values")
        if name is not None and not isinstance(name, str):
            raise TypeError("name must be a string or None")
        try:
            values = tuple(molecules)
        except TypeError as error:
            raise TypeError("molecules must be an iterable of Molecule values") from error
        if not all(isinstance(molecule, Molecule) for molecule in values):
            raise TypeError("molecules must contain only Molecule values")
        object.__setattr__(self, "_molecules", values)
        object.__setattr__(self, "_name", name or "")
        object.__setattr__(self, "_native", _chargefw._make_collection(values_native(values), name or ""))

    def __setattr__(self, name: str, value: Any) -> None:
        raise AttributeError(f"{type(self).__name__} is immutable")

    def __len__(self) -> int:
        return len(self._molecules)

    @property
    def size(self) -> int:
        return len(self._molecules)

    @property
    def empty(self) -> bool:
        return not self._molecules

    def __getitem__(self, item: int | slice) -> Molecule | tuple[Molecule, ...]:
        return self._molecules[item]

    def __iter__(self):
        return iter(self._molecules)

    @property
    def molecules(self) -> tuple[Molecule, ...]:
        return self._molecules

    @property
    def name(self) -> str:
        return self._name


def values_native(values: tuple[Molecule, ...]) -> tuple[Any, ...]:
    return tuple(value._native for value in values)


__all__ = ["__version__", "Molecule", "MoleculeCollection", "SourceIdentity"]
