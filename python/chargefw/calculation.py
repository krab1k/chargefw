"""Synchronous calculation facade."""

from __future__ import annotations

from collections.abc import Iterable, Mapping
from dataclasses import replace
from threading import Lock
from typing import cast

from ._calculation_options import RequestedCalculation
from ._calculation_values import (
    CalculationCancelledError,
    CalculationResult,
    CalculationTimings,
    ChargeFWError,
    ExecutedPlan,
    ExecutionPolicy,
    InvalidInputError,
    NoExecutablePlanError,
    NumericalFailureError,
    Plan,
    Rejection,
    _rejections,
)
from ._chargefw import calculation as _native_calculation
from ._chargefw import parameters as _native_parameters
from ._methods import Method, _method_catalog
from ._parameters import ParameterSet, ParameterSetCatalog, _parameter_set
from ._resources import default_parameter_directory
from ._types import ChargeCorrection, Execution, MethodOptionValue, ParameterMatching
from .core import Molecule, MoleculeCollection

__all__ = [
    "Assessment",
    "CalculationCancelledError",
    "CalculationResult",
    "CalculationTimings",
    "ChargeFWError",
    "ExecutedPlan",
    "ExecutionPolicy",
    "InvalidInputError",
    "NoExecutablePlanError",
    "NumericalFailureError",
    "Plan",
    "Rejection",
    "RequestedCalculation",
    "assess",
    "calculate",
    "methods",
    "parameter_sets",
]

for _value_type in (
    RequestedCalculation,
    ExecutionPolicy,
    ExecutedPlan,
    CalculationTimings,
    ChargeFWError,
    InvalidInputError,
    NoExecutablePlanError,
    NumericalFailureError,
    CalculationCancelledError,
    CalculationResult,
    Plan,
    Rejection,
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
    """Priority-ordered runnable plans and rejected alternatives for immutable molecules."""

    __slots__ = ("_native", "_molecules", "_requested", "_plans", "_rejections", "_seconds")

    def __init__(
        self,
        native: _native_calculation._NativeAssessment,
        molecules: MoleculeCollection,
        requested: RequestedCalculation,
    ) -> None:
        payload = native.report()
        self._native = native
        self._molecules = molecules
        self._requested = requested
        self._plans = tuple(
            Plan(plan, molecules, methods, parameter_sets, requested) for plan in native.plans()
        )
        self._rejections = _rejections(payload["rejections"], methods, parameter_sets)
        self._seconds = payload["applicability_seconds"]

    @property
    def plans(self) -> tuple[Plan, ...]:
        return self._plans

    @property
    def rejections(self) -> tuple[Rejection, ...]:
        return self._rejections

    @property
    def default_plan(self) -> Plan | None:
        return None if not self._plans else self._plans[0]

    @property
    def seconds(self) -> float:
        return self._seconds


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
    """Assess molecules and return reusable executable calculation plans."""

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


class _OmittedPlan:
    def __repr__(self) -> str:
        return "<automatic>"


_OMITTED_PLAN = _OmittedPlan()


def calculate(
    molecules: Molecule | MoleculeCollection | Iterable[Molecule],
    plan: Plan = cast(Plan, _OMITTED_PLAN),
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
    threads: int | None = None,
) -> CalculationResult:
    """Calculate molecules with an assessed plan, or assess and use the default plan."""

    collection = _as_collection(molecules)
    if isinstance(plan, _OmittedPlan):
        assessment = assess(
            collection,
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
            threads=0 if threads is None else threads,
        )
        if assessment.default_plan is None:
            result = CalculationResult(
                assessment._native.calculate_default(),
                collection,
                assessment._requested,
                methods,
                parameter_sets,
                assessment.rejections,
            )
            result._raise_for_status()
            return result
        return calculate(collection, assessment.default_plan, threads=threads)

    if plan is None or not isinstance(plan, Plan):
        raise TypeError("plan must be a Plan; omit it to use automatic assessment")
    if any(
        (
            method is not None,
            parameter_set is not None,
            options is not None,
            options_by_method is not None,
            parameter_matching != "strict",
            execution != "auto",
            radius is not None,
            charge_correction is not None,
            cutoff_threshold != 20_000,
            cover_threshold != 80_000,
        )
    ):
        raise TypeError("selection arguments cannot be combined with an assessed plan")
    if not plan._matches(collection):
        raise ValueError("plan was assessed for a different molecule collection")
    normalized_threads = None if threads is None else RequestedCalculation(threads=threads).threads
    requested = (
        plan._requested
        if normalized_threads is None
        else replace(plan._requested, threads=normalized_threads)
    )
    result = CalculationResult(
        plan._calculate(normalized_threads),
        collection,
        requested,
        methods,
        parameter_sets,
    )
    result._raise_for_status()
    return result


Assessment.__module__ = __name__
