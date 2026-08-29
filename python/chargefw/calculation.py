"""Synchronous calculation facade."""

from __future__ import annotations

from collections.abc import Iterable, Mapping
from threading import Lock

from ._calculation_options import RequestedCalculation
from ._calculation_values import (
    ApplicabilityReport,
    ApplicableMethod,
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
    RejectedMethod,
    _applicability_report,
    _execution_issue,
    _execution_policy,
)
from ._chargefw import calculation as _native_calculation
from ._chargefw import parameters as _native_parameters
from ._methods import ExecutionIssue, Method, _method_catalog
from ._parameters import ParameterSet, ParameterSetCatalog, _parameter_set
from ._resources import default_parameter_directory
from ._types import ChargeCorrection, Execution, MethodOptionValue, ParameterMatching
from .core import Molecule, MoleculeCollection

__all__ = [
    "ApplicableMethod",
    "ApplicabilityReport",
    "Assessment",
    "AssessmentReport",
    "CalculationCancelledError",
    "CalculationResult",
    "CalculationTimings",
    "ChargeFWError",
    "EffectiveCalculation",
    "ExecutionPolicy",
    "InvalidInputError",
    "NoExecutablePlanError",
    "NumericalFailureError",
    "RejectedMethod",
    "RequestedCalculation",
    "assess",
    "calculate",
    "methods",
    "parameter_sets",
]

for _value_type in (
    RequestedCalculation,
    ExecutionPolicy,
    ApplicableMethod,
    RejectedMethod,
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
_bundled_parameter_descriptors: ParameterSetCatalog | None = None
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


def _default_parameter_descriptors() -> ParameterSetCatalog:
    global _bundled_parameter_descriptors
    catalog = _default_parameter_catalog()
    if _bundled_parameter_descriptors is None:
        with _bundled_parameter_catalog_lock:
            if _bundled_parameter_descriptors is None:
                _bundled_parameter_descriptors = ParameterSetCatalog(
                    tuple(_parameter_set(value) for value in catalog._descriptors())
                )
    return _bundled_parameter_descriptors


parameter_sets = _default_parameter_descriptors()
methods = _method_catalog(parameter_sets)


def _as_collection(value: Molecule | MoleculeCollection | Iterable[Molecule]) -> MoleculeCollection:
    if isinstance(value, MoleculeCollection):
        return value
    if isinstance(value, Molecule):
        return MoleculeCollection((value,))
    return MoleculeCollection(value)


class Assessment:
    """One-shot applicability assessment and executable calculation plan."""

    __slots__ = ("_native", "_molecules", "_requested", "_report", "_consumed")

    _native: _native_calculation._NativeAssessment | None

    def __init__(
        self,
        native: _native_calculation._NativeAssessment,
        molecules: MoleculeCollection,
        requested: RequestedCalculation,
    ) -> None:
        payload = native.report()
        policy = payload["execution_policy"]
        self._native = native
        self._molecules = molecules
        self._requested = requested
        self._report = AssessmentReport(
            applicability=_applicability_report(payload["applicability"], methods, parameter_sets),
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
    def applicable(self) -> tuple[ApplicableMethod, ...]:
        return self._report.applicability.applicable

    @property
    def rejected(self) -> tuple[RejectedMethod, ...]:
        return self._report.applicability.rejected

    @property
    def selected(self) -> ApplicableMethod | None:
        index = self._report.applicability.selected_candidate_index
        return None if index is None else self.applicable[index]

    @property
    def execution(self) -> ExecutionPolicy | None:
        return self._report.execution_policy

    @property
    def warnings(self) -> tuple[ExecutionIssue, ...]:
        return self._report.execution_issues

    @property
    def seconds(self) -> float:
        return self._report.applicability_seconds

    @property
    def executable(self) -> bool:
        return self._report.executable

    def __enter__(self) -> Assessment:
        if self._native is None:
            raise RuntimeError("assessment is closed")
        return self

    def __exit__(self, *exception: object) -> None:
        self.close()

    def close(self) -> None:
        """Release prepared native calculation state without executing it."""

        self._native = None

    def calculate(self) -> CalculationResult:
        if self._consumed:
            raise RuntimeError("assessment can only be calculated once")
        if self._native is None:
            raise RuntimeError("assessment is closed")
        self._consumed = True
        native = self._native
        self._native = None
        result = CalculationResult(
            native.calculate(),
            self._molecules,
            self._requested,
            methods,
            parameter_sets,
        )
        result._raise_for_status()
        return result


def assess(
    molecules: Molecule | MoleculeCollection | Iterable[Molecule],
    *,
    method: str | Method | None = None,
    parameter_set: str | ParameterSet | None = None,
    options: Mapping[str, MethodOptionValue] | None = None,
    options_by_method: Mapping[str, Mapping[str, MethodOptionValue]] | None = None,
    parameter_matching: ParameterMatching = "strict",
    execution: Execution = "auto",
    radius: float | None = None,
    charge_correction: ChargeCorrection | None = None,
    cutoff_threshold: int | None = 20_000,
    cover_threshold: int | None = 80_000,
    threads: int = 0,
) -> Assessment:
    """Assess molecules and return a one-shot executable calculation plan."""

    requested = RequestedCalculation(
        method=method,
        parameter_set=parameter_set,
        options=options,
        options_by_method=options_by_method,
        parameter_matching=parameter_matching,
        execution=execution,
        radius=radius,
        charge_correction=charge_correction,
        cutoff_threshold=cutoff_threshold,
        cover_threshold=cover_threshold,
        threads=threads,
    )
    collection = _as_collection(molecules)
    native = _native_calculation._make_assessment(
        collection._native_molecules,
        collection.name,
        _default_parameter_catalog(),
        requested.method,
        requested.parameter_set,
        {
            method_id: dict(overrides)
            for method_id, overrides in requested.options_by_method.items()
        },
        requested._permissive_types,
        requested.execution,
        requested.radius,
        requested.charge_correction,
        requested.cutoff_threshold,
        requested.cover_threshold,
        requested.threads,
    )
    return Assessment(native, collection, requested)


def calculate(
    molecules: Molecule | MoleculeCollection | Iterable[Molecule],
    *,
    method: str | Method | None = None,
    parameter_set: str | ParameterSet | None = None,
    options: Mapping[str, MethodOptionValue] | None = None,
    options_by_method: Mapping[str, Mapping[str, MethodOptionValue]] | None = None,
    parameter_matching: ParameterMatching = "strict",
    execution: Execution = "auto",
    radius: float | None = None,
    charge_correction: ChargeCorrection | None = None,
    cutoff_threshold: int | None = 20_000,
    cover_threshold: int | None = 80_000,
    threads: int = 0,
) -> CalculationResult:
    """Assess and synchronously calculate charges for molecules."""

    return assess(
        molecules,
        method=method,
        parameter_set=parameter_set,
        options=options,
        options_by_method=options_by_method,
        parameter_matching=parameter_matching,
        execution=execution,
        radius=radius,
        charge_correction=charge_correction,
        cutoff_threshold=cutoff_threshold,
        cover_threshold=cover_threshold,
        threads=threads,
    ).calculate()


Assessment.__module__ = __name__
