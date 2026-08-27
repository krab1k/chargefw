"""Synchronous calculation facade and immutable result values."""

from __future__ import annotations

from collections.abc import Iterable, Mapping
from dataclasses import dataclass, field
import numbers
from operator import index as as_index
from threading import Lock
from types import MappingProxyType
from typing import Any

import numpy as np

from ._chargefw import calculation as _native_calculation
from ._chargefw import parameters as _native_parameters
from ._resources import default_parameter_directory
from .charges import ChargeAssignment
from .core import Molecule, MoleculeCollection
from .methods import (
    ExecutionAssessment,
    ExecutionIssue,
    PrerequisiteIssue,
    MethodDescriptor,
    method_descriptors,
)
from .parameters import ParameterSet, ParameterSetDescriptor, _descriptor

ExecutionSelectionKind = _native_calculation.ExecutionSelectionKind
ExecutionMode = _native_calculation.ExecutionMode
ChargeCorrectionPolicy = _native_calculation.ChargeCorrectionPolicy
ExecutionStatus = _native_calculation.ExecutionStatus
for _enum in (
    ExecutionSelectionKind,
    ExecutionMode,
    ChargeCorrectionPolicy,
    ExecutionStatus,
):
    _enum.__module__ = __name__

MethodOptionValue = bool | int | float | str

_bundled_parameter_catalog: Any | None = None
_bundled_parameter_catalog_lock = Lock()


def _default_parameter_catalog() -> Any:
    global _bundled_parameter_catalog
    with _bundled_parameter_catalog_lock:
        if _bundled_parameter_catalog is None:
            _bundled_parameter_catalog = _native_parameters._load_parameter_catalog(
                str(default_parameter_directory())
            )
        return _bundled_parameter_catalog


def _normalized_nonnegative_integer(value: Any, field_name: str) -> int:
    if isinstance(value, (bool, np.bool_)):
        raise TypeError(f"{field_name} must be an integer")
    try:
        result = as_index(value)
    except TypeError as error:
        raise TypeError(f"{field_name} must be an integer") from error
    if result < 0:
        raise ValueError(f"{field_name} must be non-negative")
    if result > np.iinfo(np.uintp).max:
        raise ValueError(f"{field_name} is outside the native size range")
    return result


def _normalized_method_options(
    values: Mapping[str, Mapping[str, MethodOptionValue]],
) -> dict[str, dict[str, MethodOptionValue]]:
    if not isinstance(values, Mapping):
        raise TypeError("method_options must be a mapping")
    result: dict[str, dict[str, MethodOptionValue]] = {}
    int_limits = np.iinfo(np.int32)
    for method_id, overrides in values.items():
        if not isinstance(method_id, str):
            raise TypeError("method option method IDs must be strings")
        if not isinstance(overrides, Mapping):
            raise TypeError(f"method_options[{method_id!r}] must be a mapping")
        normalized: dict[str, MethodOptionValue] = {}
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
                raise TypeError(
                    f"method option {method_id}.{option_id} must be bool, int, float, or str"
                )
        result[method_id] = normalized
    return result


def _frozen_method_options(
    values: Mapping[str, Mapping[str, MethodOptionValue]],
) -> Mapping[str, Mapping[str, MethodOptionValue]]:
    normalized = _normalized_method_options(values)
    return MappingProxyType(
        {
            method_id: MappingProxyType(overrides)
            for method_id, overrides in normalized.items()
        }
    )


@dataclass(frozen=True, slots=True)
class CalculationOptions:
    """Application policy for one synchronous calculation request."""

    method: str | None = None
    parameter_set_id: str | None = None
    method_options: Mapping[str, Mapping[str, MethodOptionValue]] = field(default_factory=dict)
    permissive_types: bool = False
    execution: ExecutionSelectionKind = ExecutionSelectionKind.AUTOMATIC
    radius: float | None = None
    charge_correction: ChargeCorrectionPolicy | None = None
    cutoff_atom_threshold: int | None = 20_000
    cover_atom_threshold: int | None = 80_000
    max_threads: int = 0

    def __post_init__(self) -> None:
        for field_name in ("method", "parameter_set_id"):
            value = getattr(self, field_name)
            if value is not None and not isinstance(value, str):
                raise TypeError(f"{field_name} must be a string or None")
        if not isinstance(self.permissive_types, (bool, np.bool_)):
            raise TypeError("permissive_types must be a boolean")
        if not isinstance(self.execution, ExecutionSelectionKind):
            raise TypeError("execution must be an ExecutionSelectionKind")
        if self.charge_correction is not None and not isinstance(
            self.charge_correction, ChargeCorrectionPolicy
        ):
            raise TypeError("charge_correction must be a ChargeCorrectionPolicy or None")
        if (
            self.execution is ExecutionSelectionKind.AUTOMATIC
            and self.charge_correction is not None
        ):
            raise ValueError("automatic execution does not accept a charge correction")
        if self.execution is ExecutionSelectionKind.FULL:
            if self.radius is not None:
                raise ValueError("full execution does not accept a radius")
            if self.charge_correction is not None:
                raise ValueError("full execution does not accept a charge correction")
        if self.execution in (
            ExecutionSelectionKind.CUTOFF,
            ExecutionSelectionKind.COVER,
        ) and self.radius is None:
            raise ValueError(f"{self.execution.name.lower()} execution requires a radius")
        if self.radius is not None:
            if not isinstance(self.radius, numbers.Real) or isinstance(
                self.radius, (bool, np.bool_)
            ):
                raise TypeError("radius must be a real number or None")
            if not np.isfinite(float(self.radius)) or float(self.radius) < 8.0:
                raise ValueError("radius must be finite and at least 8.0")
        for field_name in ("cutoff_atom_threshold", "cover_atom_threshold"):
            value = getattr(self, field_name)
            if value is not None:
                _normalized_nonnegative_integer(value, field_name)
        if self.cover_atom_threshold is not None and self.cutoff_atom_threshold is None:
            raise ValueError("cover_atom_threshold requires a finite cutoff_atom_threshold")
        if (
            self.cover_atom_threshold is not None
            and self.cover_atom_threshold < self.cutoff_atom_threshold
        ):
            raise ValueError("cover_atom_threshold must not be smaller than cutoff_atom_threshold")
        _normalized_nonnegative_integer(self.max_threads, "max_threads")

        object.__setattr__(self, "permissive_types", bool(self.permissive_types))
        object.__setattr__(self, "method_options", _frozen_method_options(self.method_options))
        if self.radius is not None:
            object.__setattr__(self, "radius", float(self.radius))
        if self.cutoff_atom_threshold is not None:
            object.__setattr__(self, "cutoff_atom_threshold", int(self.cutoff_atom_threshold))
        if self.cover_atom_threshold is not None:
            object.__setattr__(self, "cover_atom_threshold", int(self.cover_atom_threshold))
        object.__setattr__(self, "max_threads", int(self.max_threads))


@dataclass(frozen=True, slots=True)
class ExecutionPolicy:
    mode: ExecutionMode
    radius: float | None
    charge_correction: ChargeCorrectionPolicy


@dataclass(frozen=True, slots=True)
class ApplicableCandidate:
    method_id: str
    parameter_set_id: str | None
    execution_assessments: tuple[ExecutionAssessment, ...]


@dataclass(frozen=True, slots=True)
class RejectedCandidate:
    method_id: str
    parameter_set_id: str | None
    issues: tuple[PrerequisiteIssue, ...]


@dataclass(frozen=True, slots=True)
class ApplicabilityReport:
    applicable: tuple[ApplicableCandidate, ...]
    rejected: tuple[RejectedCandidate, ...]
    selected_candidate_index: int | None


@dataclass(frozen=True, slots=True)
class EffectiveCalculation:
    method_id: str
    parameter_set_id: str | None
    method_options: Mapping[str, MethodOptionValue]
    execution_policy: ExecutionPolicy
    execution_issues: tuple[ExecutionIssue, ...]

    def __post_init__(self) -> None:
        object.__setattr__(self, "method_options", MappingProxyType(dict(self.method_options)))


@dataclass(frozen=True, slots=True)
class CalculationTimings:
    applicability_seconds: float
    computation_seconds: float


@dataclass(frozen=True, slots=True)
class AssessmentReport:
    applicability: ApplicabilityReport
    execution_policy: ExecutionPolicy | None
    execution_issues: tuple[ExecutionIssue, ...]
    applicability_seconds: float
    executable: bool


def _prerequisite_issue(value: Mapping[str, Any]) -> PrerequisiteIssue:
    return PrerequisiteIssue(
        kind=value["kind"],
        message=value["message"],
        molecule_index=value["molecule_index"],
        atom_index=value["atom_index"],
        bond_index=value["bond_index"],
        conformer_index=value["conformer_index"],
    )


def _execution_issue(value: Mapping[str, Any]) -> ExecutionIssue:
    return ExecutionIssue(
        kind=value["kind"],
        message=value["message"],
        molecule_index=value["molecule_index"],
    )


def _execution_assessment(value: Mapping[str, Any]) -> ExecutionAssessment:
    return ExecutionAssessment(
        mode=value["mode"],
        availability=value["availability"],
        issues=tuple(_execution_issue(issue) for issue in value["issues"]),
    )


def _applicability_report(value: Mapping[str, Any]) -> ApplicabilityReport:
    return ApplicabilityReport(
        applicable=tuple(
            ApplicableCandidate(
                method_id=candidate["method_id"],
                parameter_set_id=candidate["parameter_set_id"],
                execution_assessments=tuple(
                    _execution_assessment(assessment)
                    for assessment in candidate["execution_assessments"]
                ),
            )
            for candidate in value["applicable"]
        ),
        rejected=tuple(
            RejectedCandidate(
                method_id=candidate["method_id"],
                parameter_set_id=candidate["parameter_set_id"],
                issues=tuple(_prerequisite_issue(issue) for issue in candidate["issues"]),
            )
            for candidate in value["rejected"]
        ),
        selected_candidate_index=value["selected_candidate_index"],
    )


def _execution_policy(value: Mapping[str, Any]) -> ExecutionPolicy:
    return ExecutionPolicy(
        mode=value["mode"],
        radius=value["radius"],
        charge_correction=value["charge_correction"],
    )


def _effective_calculation(value: Mapping[str, Any] | None) -> EffectiveCalculation | None:
    if value is None:
        return None
    return EffectiveCalculation(
        method_id=value["method_id"],
        parameter_set_id=value["parameter_set_id"],
        method_options=value["method_options"],
        execution_policy=_execution_policy(value["execution_policy"]),
        execution_issues=tuple(_execution_issue(issue) for issue in value["execution_issues"]),
    )


def _as_collection(value: Molecule | MoleculeCollection | Iterable[Molecule]) -> MoleculeCollection:
    if isinstance(value, MoleculeCollection):
        return value
    if isinstance(value, Molecule):
        return MoleculeCollection((value,))
    return MoleculeCollection(value)


class ChargeFWError(RuntimeError):
    """Base class for a failed ChargeFW calculation result."""

    def __init__(self, message: str, result: CalculationResult) -> None:
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
        "_status",
        "_requested",
        "_assignments",
        "_applicability",
        "_effective",
        "_failure_message",
        "_timings",
    )

    def __init__(
        self,
        payload: Mapping[str, Any],
        molecules: MoleculeCollection,
        requested: CalculationOptions,
    ) -> None:
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
                    ChargeAssignment._from_native_values(
                        item["values"],
                        molecule_index=molecule_index,
                        conformer_index=(
                            int(conformer_index) if conformer_index is not None else None
                        ),
                        source=molecule.source,
                        source_atom_ids=molecule.source_atom_ids,
                        source_conformer_id=source_conformer_id,
                    )
                )
        self._status = payload["status"]
        self._requested = requested
        self._assignments = tuple(assignments)
        self._applicability = _applicability_report(payload["applicability"])
        self._effective = _effective_calculation(payload["effective"])
        self._failure_message = payload["failure_message"]
        self._timings = CalculationTimings(
            applicability_seconds=payload["metrics"]["applicability_seconds"],
            computation_seconds=payload["metrics"]["computation_seconds"],
        )

    @property
    def status(self) -> ExecutionStatus:
        return self._status

    @property
    def assignments(self) -> tuple[ChargeAssignment, ...]:
        return self._assignments

    @property
    def applicability(self) -> ApplicabilityReport:
        return self._applicability

    @property
    def requested(self) -> CalculationOptions:
        return self._requested

    @property
    def effective(self) -> EffectiveCalculation | None:
        return self._effective

    @property
    def execution_issues(self) -> tuple[ExecutionIssue, ...]:
        if self._effective is None:
            return ()
        return self._effective.execution_issues

    @property
    def failure_message(self) -> str | None:
        return self._failure_message

    @property
    def timings(self) -> CalculationTimings:
        return self._timings

    def raise_for_status(self) -> None:
        exception_types = {
            ExecutionStatus.INVALID_INPUT_OR_REQUEST: InvalidInputError,
            ExecutionStatus.NO_EXECUTABLE_PLAN: NoExecutablePlanError,
            ExecutionStatus.NUMERICAL_FAILURE: NumericalFailureError,
            ExecutionStatus.CANCELLED: CalculationCancelledError,
        }
        exception_type = exception_types.get(self.status)
        if exception_type is not None:
            raise exception_type(self.failure_message or self.status.name, self)


class Assessment:
    """One-shot applicability assessment and executable calculation plan."""

    __slots__ = ("_native", "_molecules", "_requested", "_report", "_consumed")

    def __init__(
        self,
        native: Any,
        molecules: MoleculeCollection,
        requested: CalculationOptions,
    ) -> None:
        payload = native.report()
        policy = payload["execution_policy"]
        self._native = native
        self._molecules = molecules
        self._requested = requested
        self._report = AssessmentReport(
            applicability=_applicability_report(payload["applicability"]),
            execution_policy=_execution_policy(policy) if policy is not None else None,
            execution_issues=tuple(
                _execution_issue(issue) for issue in payload["execution_issues"]
            ),
            applicability_seconds=payload["applicability_seconds"],
            executable=payload["executable"],
        )
        self._consumed = False

    @property
    def report(self) -> AssessmentReport:
        return self._report

    @property
    def applicability(self) -> ApplicabilityReport:
        return self._report.applicability

    @property
    def execution_policy(self) -> ExecutionPolicy | None:
        return self._report.execution_policy

    @property
    def execution_issues(self) -> tuple[ExecutionIssue, ...]:
        return self._report.execution_issues

    @property
    def executable(self) -> bool:
        return self._report.executable

    def calculate(self) -> CalculationResult:
        if self._consumed:
            raise RuntimeError("assessment can only be calculated once")
        self._consumed = True
        return CalculationResult(self._native.calculate(), self._molecules, self._requested)


class Calculator:
    """Synchronous calculator using bundled or explicit immutable parameter sets."""

    __slots__ = ("_catalog", "_method_descriptors", "_parameter_set_descriptors")

    def __init__(self, parameter_sets: Iterable[ParameterSet] | None = None) -> None:
        if parameter_sets is None:
            self._catalog = _default_parameter_catalog()
        else:
            try:
                values = tuple(parameter_sets)
            except TypeError as error:
                raise TypeError(
                    "parameter_sets must be an iterable of ParameterSet values or None"
                ) from error
            if not values:
                raise ValueError("parameter_sets must not be empty")
            if not all(isinstance(value, ParameterSet) for value in values):
                raise TypeError("parameter_sets must contain only ParameterSet values")
            ids = [value.id for value in values]
            if len(ids) != len(set(ids)):
                raise ValueError("parameter_sets must have unique IDs")
            self._catalog = _native_parameters._load_parameter_catalog_from_sets(
                [value._native for value in values]
            )
        self._method_descriptors = method_descriptors()
        self._parameter_set_descriptors = tuple(
            _descriptor(value) for value in self._catalog._descriptors()
        )

    @property
    def methods(self) -> tuple[MethodDescriptor, ...]:
        """Value-only descriptors for the built-in method registry."""

        return self._method_descriptors

    @property
    def parameter_sets(self) -> tuple[ParameterSetDescriptor, ...]:
        """Value-only descriptors for this calculator's parameter catalog."""

        return self._parameter_set_descriptors

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
        native = _native_calculation._make_assessment(
            collection._native,
            self._catalog,
            options.method,
            options.parameter_set_id,
            {
                method_id: dict(overrides)
                for method_id, overrides in options.method_options.items()
            },
            options.permissive_types,
            options.execution,
            options.radius,
            options.charge_correction,
            options.cutoff_atom_threshold,
            options.cover_atom_threshold,
            options.max_threads,
        )
        return Assessment(native, collection, options)

    def calculate(
        self,
        molecules: Molecule | MoleculeCollection | Iterable[Molecule],
        options: CalculationOptions | None = None,
    ) -> CalculationResult:
        return self.assess(molecules, options).calculate()
