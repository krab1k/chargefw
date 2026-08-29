"""Private immutable calculation result values and native-payload conversion."""

from __future__ import annotations

from collections.abc import Mapping
from dataclasses import dataclass
from types import MappingProxyType
from typing import cast

from ._calculation_options import RequestedCalculation
from ._methods import (
    ExecutionIssue,
    ExecutionSupport,
    Method,
    MethodCatalog,
    PrerequisiteIssue,
)
from ._parameters import ParameterSet, ParameterSetCatalog
from ._payloads import (
    ApplicabilityReportPayload,
    EffectiveCalculationPayload,
    ExecutionAssessmentPayload,
    ExecutionIssuePayload,
    ExecutionPolicyPayload,
    ExecutionResultPayload,
    PrerequisiteIssuePayload,
)
from ._types import (
    ChargeCorrection,
    ExecutionAvailability,
    ExecutionIssueKind,
    ExecutionMode,
    ExecutionStatus,
    MethodOptionValue,
    PrerequisiteIssueKind,
)
from .charges import ChargeAssignment
from .core import MoleculeCollection


@dataclass(frozen=True, slots=True)
class ExecutionPolicy:
    mode: ExecutionMode
    radius: float | None
    charge_correction: ChargeCorrection


@dataclass(frozen=True, slots=True)
class ApplicableMethod:
    method: Method
    parameter_set: ParameterSet | None
    executions: tuple[ExecutionSupport, ...]


@dataclass(frozen=True, slots=True)
class RejectedMethod:
    method: Method
    parameter_set: ParameterSet | None
    issues: tuple[PrerequisiteIssue, ...]


@dataclass(frozen=True, slots=True)
class ApplicabilityReport:
    applicable: tuple[ApplicableMethod, ...]
    rejected: tuple[RejectedMethod, ...]
    selected_candidate_index: int | None


@dataclass(frozen=True, slots=True)
class EffectiveCalculation:
    method: Method
    parameter_set: ParameterSet | None
    options: Mapping[str, MethodOptionValue]
    execution: ExecutionPolicy
    warnings: tuple[ExecutionIssue, ...]

    def __post_init__(self) -> None:
        object.__setattr__(self, "options", MappingProxyType(dict(self.options)))


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


def _prerequisite_issue(value: PrerequisiteIssuePayload) -> PrerequisiteIssue:
    return PrerequisiteIssue(
        kind=cast(PrerequisiteIssueKind, value["kind"].name.lower()),
        message=value["message"],
        molecule_index=value["molecule_index"],
        atom_index=value["atom_index"],
        bond_index=value["bond_index"],
        conformer_index=value["conformer_index"],
    )


def _execution_issue(value: ExecutionIssuePayload) -> ExecutionIssue:
    return ExecutionIssue(
        kind=cast(ExecutionIssueKind, value["kind"].name.lower()),
        message=value["message"],
        molecule_index=value["molecule_index"],
    )


def _execution_support(value: ExecutionAssessmentPayload) -> ExecutionSupport:
    return ExecutionSupport(
        mode=cast(ExecutionMode, value["mode"].name.lower()),
        availability=cast(ExecutionAvailability, value["availability"].name.lower()),
        issues=tuple(_execution_issue(issue) for issue in value["issues"]),
    )


def _resolved_parameter_set(
    id: str | None, parameter_sets: ParameterSetCatalog
) -> ParameterSet | None:
    return None if id is None else parameter_sets[id]


def _applicability_report(
    value: ApplicabilityReportPayload,
    methods: MethodCatalog,
    parameter_sets: ParameterSetCatalog,
) -> ApplicabilityReport:
    return ApplicabilityReport(
        applicable=tuple(
            ApplicableMethod(
                method=methods[candidate["method_id"]],
                parameter_set=_resolved_parameter_set(
                    candidate["parameter_set_id"], parameter_sets
                ),
                executions=tuple(
                    _execution_support(assessment)
                    for assessment in candidate["execution_assessments"]
                ),
            )
            for candidate in value["applicable"]
        ),
        rejected=tuple(
            RejectedMethod(
                method=methods[candidate["method_id"]],
                parameter_set=_resolved_parameter_set(
                    candidate["parameter_set_id"], parameter_sets
                ),
                issues=tuple(_prerequisite_issue(issue) for issue in candidate["issues"]),
            )
            for candidate in value["rejected"]
        ),
        selected_candidate_index=value["selected_candidate_index"],
    )


def _execution_policy(value: ExecutionPolicyPayload) -> ExecutionPolicy:
    return ExecutionPolicy(
        mode=cast(ExecutionMode, value["mode"].name.lower()),
        radius=value["radius"],
        charge_correction=cast(ChargeCorrection, value["charge_correction"].name.lower()),
    )


def _effective_calculation(
    value: EffectiveCalculationPayload | None,
    methods: MethodCatalog,
    parameter_sets: ParameterSetCatalog,
) -> EffectiveCalculation | None:
    if value is None:
        return None
    return EffectiveCalculation(
        method=methods[value["method_id"]],
        parameter_set=_resolved_parameter_set(value["parameter_set_id"], parameter_sets),
        options=value["method_options"],
        execution=_execution_policy(value["execution_policy"]),
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
        "_applicability",
        "_effective",
        "_failure_message",
        "_timings",
    )

    _status: ExecutionStatus
    _requested: RequestedCalculation
    _assignments: tuple[ChargeAssignment, ...]
    _applicability: ApplicabilityReport
    _effective: EffectiveCalculation | None
    _failure_message: str | None
    _timings: CalculationTimings

    def __init__(
        self,
        payload: ExecutionResultPayload,
        molecules: MoleculeCollection,
        requested: RequestedCalculation,
        methods: MethodCatalog,
        parameter_sets: ParameterSetCatalog,
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
        self._status = cast(ExecutionStatus, payload["status"].name.lower())
        self._requested = requested
        self._assignments = tuple(assignments)
        self._applicability = _applicability_report(
            payload["applicability"], methods, parameter_sets
        )
        self._effective = _effective_calculation(payload["effective"], methods, parameter_sets)
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
    def applicable(self) -> tuple[ApplicableMethod, ...]:
        return self._applicability.applicable

    @property
    def rejected(self) -> tuple[RejectedMethod, ...]:
        return self._applicability.rejected

    @property
    def selected(self) -> ApplicableMethod | None:
        index = self._applicability.selected_candidate_index
        return None if index is None else self.applicable[index]

    @property
    def requested(self) -> RequestedCalculation:
        return self._requested

    @property
    def effective(self) -> EffectiveCalculation | None:
        return self._effective

    @property
    def warnings(self) -> tuple[ExecutionIssue, ...]:
        if self._effective is None:
            return ()
        return self._effective.warnings

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
