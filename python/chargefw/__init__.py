"""ChargeFW Python API."""

from __future__ import annotations

from collections.abc import Hashable, Iterable, Mapping, Sequence
from copy import deepcopy
from dataclasses import dataclass, field
import importlib.resources
import numbers
from operator import index as as_index
from pathlib import Path
from typing import Any
from types import MappingProxyType

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


def _normalized_nonnegative_integer(value: Any, field: str) -> int:
    try:
        result = as_index(value)
    except TypeError as error:
        raise TypeError(f"{field} must be an integer") from error
    if result < 0:
        raise ValueError(f"{field} must be non-negative")
    return result


def _normalized_method_options(
    values: Mapping[str, Mapping[str, bool | int | float | str]]
) -> dict[str, dict[str, bool | int | float | str]]:
    if not isinstance(values, Mapping):
        raise TypeError("method_options must be a mapping")
    result: dict[str, dict[str, bool | int | float | str]] = {}
    int_limits = np.iinfo(np.int32)
    for method_id, overrides in values.items():
        if not isinstance(method_id, str):
            raise TypeError("method option method IDs must be strings")
        if not isinstance(overrides, Mapping):
            raise TypeError(f"method_options[{method_id!r}] must be a mapping")
        normalized: dict[str, bool | int | float | str] = {}
        for option_id, value in overrides.items():
            if not isinstance(option_id, str):
                raise TypeError("method option IDs must be strings")
            if isinstance(value, (bool, np.bool_)):
                normalized[option_id] = bool(value)
            elif isinstance(value, numbers.Integral):
                integer = int(value)
                if integer < int_limits.min or integer > int_limits.max:
                    raise ValueError(f"method option {method_id}.{option_id} is outside int range")
                normalized[option_id] = integer
            elif isinstance(value, numbers.Real):
                floating = float(value)
                if not np.isfinite(floating):
                    raise ValueError(f"method option {method_id}.{option_id} must be finite")
                normalized[option_id] = floating
            elif isinstance(value, str):
                normalized[option_id] = value
            else:
                raise TypeError(f"method option {method_id}.{option_id} must be bool, int, float, or str")
        result[method_id] = normalized
    return result


@dataclass(frozen=True, slots=True)
class CalculationOptions:
    """Application policy for one synchronous calculation request."""

    method: str | None = None
    parameter_set_id: str | None = None
    method_options: Mapping[str, Mapping[str, bool | int | float | str]] = field(
        default_factory=dict
    )
    permissive_types: bool = False
    execution: str = "auto"
    radius: float | None = None
    charge_correction: str | None = None
    cutoff_atom_threshold: int | None = 20_000
    cover_atom_threshold: int | None = 80_000
    max_threads: int = 0

    def __post_init__(self) -> None:
        for field in ("method", "parameter_set_id"):
            value = getattr(self, field)
            if value is not None and not isinstance(value, str):
                raise TypeError(f"{field} must be a string or None")
        if not isinstance(self.permissive_types, (bool, np.bool_)):
            raise TypeError("permissive_types must be a boolean")
        if self.execution not in {"auto", "full", "cutoff", "cover"}:
            raise ValueError("execution must be 'auto', 'full', 'cutoff', or 'cover'")
        if self.charge_correction not in {None, "none", "uniform"}:
            raise ValueError("charge_correction must be None, 'none', or 'uniform'")
        if self.execution == "auto" and self.charge_correction is not None:
            raise ValueError("automatic execution does not accept a charge correction")
        if self.execution == "full":
            if self.radius is not None:
                raise ValueError("full execution does not accept a radius")
            if self.charge_correction is not None:
                raise ValueError("full execution does not accept a charge correction")
        if self.execution in {"cutoff", "cover"} and self.radius is None:
            raise ValueError(f"{self.execution} execution requires a radius")
        if self.radius is not None:
            if not isinstance(self.radius, numbers.Real) or isinstance(self.radius, (bool, np.bool_)):
                raise TypeError("radius must be a real number or None")
            if not np.isfinite(float(self.radius)) or float(self.radius) < 8.0:
                raise ValueError("radius must be finite and at least 8.0")
        for field in ("cutoff_atom_threshold", "cover_atom_threshold"):
            value = getattr(self, field)
            if value is not None:
                _normalized_nonnegative_integer(value, field)
        if self.cover_atom_threshold is not None and self.cutoff_atom_threshold is None:
            raise ValueError("cover_atom_threshold requires a finite cutoff_atom_threshold")
        if (
            self.cover_atom_threshold is not None
            and self.cover_atom_threshold < self.cutoff_atom_threshold
        ):
            raise ValueError("cover_atom_threshold must not be smaller than cutoff_atom_threshold")
        _normalized_nonnegative_integer(self.max_threads, "max_threads")
        object.__setattr__(self, "permissive_types", bool(self.permissive_types))
        normalized_options = _normalized_method_options(self.method_options)
        object.__setattr__(
            self,
            "method_options",
            MappingProxyType(
                {method_id: MappingProxyType(overrides) for method_id, overrides in normalized_options.items()}
            ),
        )
        if self.radius is not None:
            object.__setattr__(self, "radius", float(self.radius))
        if self.cutoff_atom_threshold is not None:
            object.__setattr__(
                self, "cutoff_atom_threshold", int(self.cutoff_atom_threshold)
            )
        if self.cover_atom_threshold is not None:
            object.__setattr__(self, "cover_atom_threshold", int(self.cover_atom_threshold))
        object.__setattr__(self, "max_threads", int(self.max_threads))


def _default_parameter_directory() -> str:
    package_resources = importlib.resources.files("chargefw").joinpath("_data", "parameters")
    if package_resources.is_dir():
        return str(package_resources)
    return str(Path(__file__).resolve().parents[2] / "data" / "parameters")


def _as_collection(value: Molecule | MoleculeCollection | Iterable[Molecule]) -> MoleculeCollection:
    if isinstance(value, MoleculeCollection):
        return value
    if isinstance(value, Molecule):
        return MoleculeCollection((value,))
    return MoleculeCollection(value)


def _requested_options(options: CalculationOptions) -> dict[str, Any]:
    return {
        "method_id": options.method,
        "parameter_set_id": options.parameter_set_id,
        "method_options": {
            method_id: dict(overrides) for method_id, overrides in options.method_options.items()
        },
        "permissive_types": options.permissive_types,
        "execution": options.execution,
        "radius": options.radius,
        "charge_correction": options.charge_correction,
        "cutoff_atom_threshold": options.cutoff_atom_threshold,
        "cover_atom_threshold": options.cover_atom_threshold,
        "max_threads": options.max_threads,
    }


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


class ChargeFWError(RuntimeError):
    """Base class for a failed ChargeFW calculation result."""

    def __init__(self, message: str, result: "CalculationResult") -> None:
        super().__init__(message)
        self.result = result


class InvalidInputError(ChargeFWError):
    pass


class NoExecutablePlanError(ChargeFWError):
    pass


class NumericalFailureError(ChargeFWError):
    pass


class CalculationCancelledError(ChargeFWError):
    pass


class CalculationResult:
    """Owned calculation status, assignments, diagnostics, and provenance."""

    __slots__ = (
        "_payload",
        "_molecules",
        "_requested",
        "_assignments",
        "_applicability",
        "_effective",
    )

    def __init__(
        self,
        payload: Mapping[str, Any],
        molecules: MoleculeCollection,
        requested: Mapping[str, Any],
    ) -> None:
        self._payload = deepcopy(dict(payload))
        self._molecules = molecules
        self._requested = deepcopy(dict(requested))
        assignments: list[ChargeAssignment] = []
        charges = payload.get("charges")
        if charges is not None:
            for item in charges["assignments"]:
                molecule_index = int(item["molecule_index"])
                conformer_index = item["conformer_index"]
                molecule = molecules[molecule_index]
                source_conformer_id = (
                    molecule.conformer_ids[int(conformer_index)]
                    if conformer_index is not None
                    else None
                )
                assignments.append(
                    ChargeAssignment(
                        values=np.asarray(item["values"], dtype=np.float64),
                        molecule_index=molecule_index,
                        conformer_index=(
                            int(conformer_index) if conformer_index is not None else None
                        ),
                        source=molecule.source,
                        source_atom_ids=molecule.source_atom_ids,
                        source_conformer_id=source_conformer_id,
                    )
                )
        self._assignments = tuple(assignments)
        self._applicability = deepcopy(payload["applicability"])
        self._effective = deepcopy(payload["effective"])

    @property
    def status(self) -> str:
        return self._payload["status"]

    @property
    def assignments(self) -> tuple[ChargeAssignment, ...]:
        return self._assignments

    @property
    def applicability(self) -> dict[str, Any]:
        return deepcopy(self._applicability)

    @property
    def requested(self) -> dict[str, Any]:
        return deepcopy(self._requested)

    @property
    def effective(self) -> dict[str, Any] | None:
        return deepcopy(self._effective)

    @property
    def execution_issues(self) -> list[dict[str, Any]]:
        if self._effective is None:
            return []
        return deepcopy(self._effective["execution_issues"])

    @property
    def failure_message(self) -> str | None:
        return self._payload["failure_message"]

    @property
    def timings(self) -> dict[str, float]:
        return deepcopy(self._payload["metrics"])

    def raise_for_status(self) -> None:
        exception_types = {
            "invalid_input_or_request": InvalidInputError,
            "no_executable_plan": NoExecutablePlanError,
            "numerical_failure": NumericalFailureError,
            "cancelled": CalculationCancelledError,
        }
        exception_type = exception_types.get(self.status)
        if exception_type is not None:
            raise exception_type(self.failure_message or self.status, self)


class Assessment:
    """One-shot applicability assessment and executable calculation plan."""

    __slots__ = ("_native", "_molecules", "_requested", "_report", "_consumed")

    def __init__(
        self,
        native: Any,
        molecules: MoleculeCollection,
        requested: Mapping[str, Any],
    ) -> None:
        self._native = native
        self._molecules = molecules
        self._requested = deepcopy(dict(requested))
        self._report = deepcopy(native.report())
        self._consumed = False

    @property
    def report(self) -> dict[str, Any]:
        return deepcopy(self._report)

    @property
    def applicability(self) -> dict[str, Any]:
        return deepcopy(self._report["applicability"])

    @property
    def execution_policy(self) -> dict[str, Any] | None:
        return deepcopy(self._report["execution_policy"])

    @property
    def execution_issues(self) -> list[dict[str, Any]]:
        return deepcopy(self._report["execution_issues"])

    @property
    def executable(self) -> bool:
        return bool(self._report["executable"])

    def calculate(self) -> CalculationResult:
        if self._consumed:
            raise RuntimeError("assessment can only be calculated once")
        self._consumed = True
        return CalculationResult(self._native.calculate(), self._molecules, self._requested)


class Calculator:
    """Synchronous calculator using the bundled immutable parameter catalog."""

    __slots__ = ("_catalog",)

    def __init__(self, parameter_directory: str | Path | None = None) -> None:
        directory = _default_parameter_directory() if parameter_directory is None else parameter_directory
        if not isinstance(directory, (str, Path)):
            raise TypeError("parameter_directory must be a path or None")
        self._catalog = _chargefw._load_parameter_catalog(str(directory))

    def assess(
        self,
        molecules: Molecule | MoleculeCollection | Iterable[Molecule],
        options: CalculationOptions | None = None,
    ) -> Assessment:
        if options is None:
            options = CalculationOptions()
        if not isinstance(options, CalculationOptions):
            raise TypeError("options must be a CalculationOptions value or None")
        collection = _as_collection(molecules)
        native = _chargefw._make_assessment(
            collection._native,
            self._catalog,
            options.method,
            options.parameter_set_id,
            dict(options.method_options),
            options.permissive_types,
            options.execution,
            options.radius,
            options.charge_correction,
            options.cutoff_atom_threshold,
            options.cover_atom_threshold,
            options.max_threads,
        )
        return Assessment(native, collection, _requested_options(options))

    def calculate(
        self,
        molecules: Molecule | MoleculeCollection | Iterable[Molecule],
        options: CalculationOptions | None = None,
    ) -> CalculationResult:
        return self.assess(molecules, options).calculate()


__all__ = [
    "__version__",
    "Molecule",
    "MoleculeCollection",
    "SourceIdentity",
    "CalculationOptions",
    "ChargeAssignment",
    "CalculationResult",
    "Assessment",
    "Calculator",
    "ChargeFWError",
    "InvalidInputError",
    "NoExecutablePlanError",
    "NumericalFailureError",
    "CalculationCancelledError",
]
