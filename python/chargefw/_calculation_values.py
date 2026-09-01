"""Private immutable calculation result values and native-payload conversion."""

from __future__ import annotations

from collections.abc import Mapping
from dataclasses import dataclass
from types import MappingProxyType

from ._calculation_options import RequestedCalculation
from ._chargefw import calculation as _native_calculation
from ._methods import (
    ExecutionIssue,
    Method,
    MethodCatalog,
    PrerequisiteIssue,
)
from ._parameters import ParameterSet, ParameterSetCatalog
from ._payloads import (
    EffectiveCalculationPayload,
    ExecutionIssuePayload,
    ExecutionPolicyPayload,
    ExecutionResultPayload,
    PrerequisiteIssuePayload,
    RejectionPayload,
)
from ._types import (
    ChargeCorrection,
    ExecutionMode,
    ExecutionStatus,
    MethodOptionValue,
)
from .charges import ChargeAssignment
from .core import MoleculeCollection


@dataclass(frozen=True, slots=True)
class ExecutionPolicy:
    mode: ExecutionMode
    radius: float | None
    charge_correction: ChargeCorrection


class Plan:
    """A complete reusable execution plan bound to one assessed molecule collection."""

    __slots__ = (
        "_method",
        "_parameter_set",
        "_options",
        "_policy",
        "_warnings",
        "_native",
        "_molecules",
        "_requested",
    )

    def __init__(
        self,
        native: _native_calculation._NativePlan,
        molecules: MoleculeCollection,
        methods: MethodCatalog,
        parameter_sets: ParameterSetCatalog,
        requested: RequestedCalculation,
    ) -> None:
        payload = native.report()
        self._method = methods[payload["method_id"]]
        self._parameter_set = _resolved_parameter_set(payload["parameter_set_id"], parameter_sets)
        self._options = MappingProxyType(dict(payload["method_options"]))
        self._policy = _execution_policy(payload["execution_policy"])
        self._warnings = tuple(_execution_issue(issue) for issue in payload["warnings"])
        self._native = native
        self._molecules = molecules
        self._requested = requested

    @property
    def method(self) -> Method:
        return self._method

    @property
    def parameter_set(self) -> ParameterSet | None:
        return self._parameter_set

    @property
    def options(self) -> Mapping[str, MethodOptionValue]:
        return self._options

    @property
    def policy(self) -> ExecutionPolicy:
        return self._policy

    @property
    def warnings(self) -> tuple[ExecutionIssue, ...]:
        return self._warnings

    def _matches(self, molecules: MoleculeCollection) -> bool:
        return (
            molecules.name == self._molecules.name
            and len(molecules) == len(self._molecules)
            and all(
                supplied is assessed
                for supplied, assessed in zip(molecules, self._molecules, strict=True)
            )
        )

    def _calculate(self, threads: int | None) -> ExecutionResultPayload:
        return self._native.calculate(threads)

    def __repr__(self) -> str:
        parameter_set = self.parameter_set.id if self.parameter_set is not None else None
        return (
            f"{type(self).__name__}(method={self.method.id!r}, "
            f"parameter_set={parameter_set!r}, mode={self.policy.mode!r})"
        )


@dataclass(frozen=True, slots=True)
class Rejection:
    method: Method
    parameter_set: ParameterSet | None
    policy: ExecutionPolicy | None
    issues: tuple[PrerequisiteIssue | ExecutionIssue, ...]


@dataclass(frozen=True, slots=True)
class ExecutedPlan:
    method: Method
    parameter_set: ParameterSet | None
    options: Mapping[str, MethodOptionValue]
    policy: ExecutionPolicy
    warnings: tuple[ExecutionIssue, ...]

    def __post_init__(self) -> None:
        object.__setattr__(self, "options", MappingProxyType(dict(self.options)))


@dataclass(frozen=True, slots=True)
class CalculationTimings:
    applicability_seconds: float
    computation_seconds: float


def _prerequisite_issue(value: PrerequisiteIssuePayload) -> PrerequisiteIssue:
    return PrerequisiteIssue(
        kind=value["kind"],
        message=value["message"],
        molecule_index=value["molecule_index"],
        atom_index=value["atom_index"],
        bond_index=value["bond_index"],
        conformer_index=value["conformer_index"],
    )


def _execution_issue(value: ExecutionIssuePayload) -> ExecutionIssue:
    return ExecutionIssue(
        kind=value["kind"],
        message=value["message"],
        molecule_index=value["molecule_index"],
    )


def _resolved_parameter_set(
    id: str | None, parameter_sets: ParameterSetCatalog
) -> ParameterSet | None:
    return None if id is None else parameter_sets[id]


def _rejections(
    values: list[RejectionPayload],
    methods: MethodCatalog,
    parameter_sets: ParameterSetCatalog,
) -> tuple[Rejection, ...]:
    return tuple(
        Rejection(
            method=methods[value["method_id"]],
            parameter_set=_resolved_parameter_set(value["parameter_set_id"], parameter_sets),
            policy=(
                None
                if value["execution_policy"] is None
                else _execution_policy(value["execution_policy"])
            ),
            issues=(
                *(_prerequisite_issue(issue) for issue in value["prerequisite_issues"]),
                *(_execution_issue(issue) for issue in value["execution_issues"]),
            ),
        )
        for value in values
    )


def _execution_policy(value: ExecutionPolicyPayload) -> ExecutionPolicy:
    return ExecutionPolicy(
        mode=value["mode"],
        radius=value["radius"],
        charge_correction=value["charge_correction"],
    )


def _executed_plan(
    value: EffectiveCalculationPayload | None,
    methods: MethodCatalog,
    parameter_sets: ParameterSetCatalog,
) -> ExecutedPlan | None:
    if value is None:
        return None
    return ExecutedPlan(
        method=methods[value["method_id"]],
        parameter_set=_resolved_parameter_set(value["parameter_set_id"], parameter_sets),
        options=value["method_options"],
        policy=_execution_policy(value["execution_policy"]),
        warnings=tuple(_execution_issue(issue) for issue in value["execution_issues"]),
    )


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
        "_rejections",
        "_plan",
        "_failure_message",
        "_timings",
    )

    _status: ExecutionStatus
    _requested: RequestedCalculation
    _assignments: tuple[ChargeAssignment, ...]
    _rejections: tuple[Rejection, ...]
    _plan: ExecutedPlan | None
    _failure_message: str | None
    _timings: CalculationTimings

    def __init__(
        self,
        payload: ExecutionResultPayload,
        molecules: MoleculeCollection,
        requested: RequestedCalculation,
        methods: MethodCatalog,
        parameter_sets: ParameterSetCatalog,
        rejections: tuple[Rejection, ...] = (),
    ) -> None:
        assignments: list[ChargeAssignment] = []
        charges = payload["charges"]
        if charges is not None:
            for item in charges["assignments"]:
                molecule_index = int(item["molecule_index"])
                conformer_index = item["conformer_index"]
                molecule = molecules[molecule_index]
                conformer_id = (
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
                        atom_ids=molecule.atom_ids,
                        conformer_id=conformer_id,
                    )
                )
        self._status = payload["status"]
        self._requested = requested
        self._assignments = tuple(assignments)
        self._rejections = rejections or _rejections(payload["rejections"], methods, parameter_sets)
        self._plan = _executed_plan(payload["effective"], methods, parameter_sets)
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
    def rejections(self) -> tuple[Rejection, ...]:
        return self._rejections

    @property
    def requested(self) -> RequestedCalculation:
        return self._requested

    @property
    def plan(self) -> ExecutedPlan | None:
        return self._plan

    @property
    def warnings(self) -> tuple[ExecutionIssue, ...]:
        if self._plan is None:
            return ()
        return self._plan.warnings

    @property
    def failure_message(self) -> str | None:
        return self._failure_message

    @property
    def timings(self) -> CalculationTimings:
        return self._timings

    def _raise_for_status(self) -> None:
        exception_types = {
            "invalid_input_or_request": InvalidInputError,
            "no_executable_plan": NoExecutablePlanError,
            "numerical_failure": NumericalFailureError,
            "cancelled": CalculationCancelledError,
        }
        exception_type = exception_types.get(self.status)
        if exception_type is not None:
            raise exception_type(self.failure_message or self.status, self)
