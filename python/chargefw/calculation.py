"""Synchronous calculation facade."""

from __future__ import annotations

from collections.abc import Iterable
from threading import Lock

from ._calculation_options import CalculationOptions
from ._calculation_values import (
    ApplicabilityReport,
    ApplicableCandidate,
    AssessmentReport,
    CalculationCancelledError,
    CalculationResult,
    CalculationTimings,
    ChargeFWError,
    EffectiveCalculation,
    ExecutionPolicy,
    InvalidInputError,
    NoExecutablePlanError,
    NumericalFailureError,
    RejectedCandidate,
    _applicability_report,
    _execution_issue,
    _execution_policy,
)
from ._chargefw import calculation as _native_calculation
from ._chargefw import parameters as _native_parameters
from ._resources import default_parameter_directory
from .core import Molecule, MoleculeCollection
from .methods import ExecutionIssue, MethodDescriptor, method_descriptors
from .parameters import ParameterSet, ParameterSetDescriptor, _descriptor

__all__ = [
    "ApplicableCandidate",
    "ApplicabilityReport",
    "Assessment",
    "AssessmentReport",
    "CalculationCancelledError",
    "CalculationOptions",
    "CalculationResult",
    "CalculationTimings",
    "Calculator",
    "ChargeCorrectionPolicy",
    "ChargeFWError",
    "EffectiveCalculation",
    "ExecutionMode",
    "ExecutionPolicy",
    "ExecutionSelectionKind",
    "ExecutionStatus",
    "InvalidInputError",
    "NoExecutablePlanError",
    "NumericalFailureError",
    "RejectedCandidate",
]

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

for _value_type in (
    CalculationOptions,
    ExecutionPolicy,
    ApplicableCandidate,
    RejectedCandidate,
    ApplicabilityReport,
    EffectiveCalculation,
    CalculationTimings,
    AssessmentReport,
    ChargeFWError,
    InvalidInputError,
    NoExecutablePlanError,
    NumericalFailureError,
    CalculationCancelledError,
    CalculationResult,
):
    _value_type.__module__ = __name__

_bundled_parameter_catalog: _native_parameters._NativeParameterCatalog | None = None
_bundled_parameter_descriptors: tuple[ParameterSetDescriptor, ...] | None = None
_bundled_parameter_catalog_lock = Lock()


def _default_parameter_catalog() -> _native_parameters._NativeParameterCatalog:
    global _bundled_parameter_catalog
    if _bundled_parameter_catalog is None:
        with _bundled_parameter_catalog_lock:
            if _bundled_parameter_catalog is not None:
                return _bundled_parameter_catalog
            _bundled_parameter_catalog = _native_parameters._load_parameter_catalog(
                str(default_parameter_directory())
            )
    return _bundled_parameter_catalog


def _default_parameter_descriptors() -> tuple[ParameterSetDescriptor, ...]:
    global _bundled_parameter_descriptors
    catalog = _default_parameter_catalog()
    if _bundled_parameter_descriptors is None:
        with _bundled_parameter_catalog_lock:
            if _bundled_parameter_descriptors is None:
                _bundled_parameter_descriptors = tuple(
                    _descriptor(value) for value in catalog._descriptors()
                )
    return _bundled_parameter_descriptors


def _as_collection(value: Molecule | MoleculeCollection | Iterable[Molecule]) -> MoleculeCollection:
    if isinstance(value, MoleculeCollection):
        return value
    if isinstance(value, Molecule):
        return MoleculeCollection((value,))
    return MoleculeCollection(value)


class Assessment:
    """One-shot applicability assessment and executable calculation plan."""

    __slots__ = ("_native", "_molecules", "_requested", "_report", "_consumed")

    def __init__(
        self,
        native: _native_calculation._NativeAssessment,
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
            parameter_set_descriptors = _default_parameter_descriptors()
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
            parameter_set_descriptors = tuple(
                _descriptor(value) for value in self._catalog._descriptors()
            )
        self._method_descriptors = method_descriptors()
        self._parameter_set_descriptors = parameter_set_descriptors

    @property
    def methods(self) -> tuple[MethodDescriptor, ...]:
        """Value-only descriptors for the built-in method registry."""

        return self._method_descriptors

    @property
    def parameter_sets(self) -> tuple[ParameterSetDescriptor, ...]:
        """Value-only descriptors for this calculator's parameter catalog."""

        return self._parameter_set_descriptors

    def __repr__(self) -> str:
        return f"{type(self).__name__}(parameter_sets={len(self.parameter_sets)})"

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
            collection._native_molecules,
            collection.name,
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


Assessment.__module__ = __name__
Calculator.__module__ = __name__
