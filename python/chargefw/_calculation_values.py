"""Private immutable calculation result values and native-payload conversion."""

from __future__ import annotations

from collections.abc import Mapping
from dataclasses import dataclass
from types import MappingProxyType

from ._calculation_options import CalculationOptions
from ._chargefw import calculation as _native_calculation
from ._payloads import (
    ApplicabilityReportPayload,
    EffectiveCalculationPayload,
    ExecutionAssessmentPayload,
    ExecutionIssuePayload,
    ExecutionPolicyPayload,
    ExecutionResultPayload,
    PrerequisiteIssuePayload,
)
from ._types import MethodOptionValue
from .charges import ChargeAssignment
from .core import MoleculeCollection
from ._methods import ExecutionAssessment, ExecutionIssue, PrerequisiteIssue


@dataclass(frozen=True, slots=True)
class ExecutionPolicy:
    mode: _native_calculation.ExecutionMode
    radius: float | None
    charge_correction: _native_calculation.ChargeCorrectionPolicy


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


def _execution_assessment(value: ExecutionAssessmentPayload) -> ExecutionAssessment:
    return ExecutionAssessment(
        mode=value["mode"],
        availability=value["availability"],
        issues=tuple(_execution_issue(issue) for issue in value["issues"]),
    )


def _applicability_report(value: ApplicabilityReportPayload) -> ApplicabilityReport:
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


def _execution_policy(value: ExecutionPolicyPayload) -> ExecutionPolicy:
    return ExecutionPolicy(
        mode=value["mode"],
        radius=value["radius"],
        charge_correction=value["charge_correction"],
    )


def _effective_calculation(value: EffectiveCalculationPayload | None) -> EffectiveCalculation | None:
    if value is None:
        return None
    return EffectiveCalculation(
        method_id=value["method_id"],
        parameter_set_id=value["parameter_set_id"],
        method_options=value["method_options"],
        execution_policy=_execution_policy(value["execution_policy"]),
        execution_issues=tuple(_execution_issue(issue) for issue in value["execution_issues"]),
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

    _status: _native_calculation.ExecutionStatus
    _requested: CalculationOptions
    _assignments: tuple[ChargeAssignment, ...]
    _applicability: ApplicabilityReport
    _effective: EffectiveCalculation | None
    _failure_message: str | None
    _timings: CalculationTimings

    def __init__(
        self,
        payload: ExecutionResultPayload,
        molecules: MoleculeCollection,
        requested: CalculationOptions,
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
        self._applicability = _applicability_report(payload["applicability"])
        self._effective = _effective_calculation(payload["effective"])
        self._failure_message = payload["failure_message"]
        self._timings = CalculationTimings(
            applicability_seconds=payload["metrics"]["applicability_seconds"],
            computation_seconds=payload["metrics"]["computation_seconds"],
        )

    @property
    def status(self) -> _native_calculation.ExecutionStatus:
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
            _native_calculation.ExecutionStatus.INVALID_INPUT_OR_REQUEST: InvalidInputError,
            _native_calculation.ExecutionStatus.NO_EXECUTABLE_PLAN: NoExecutablePlanError,
            _native_calculation.ExecutionStatus.NUMERICAL_FAILURE: NumericalFailureError,
            _native_calculation.ExecutionStatus.CANCELLED: CalculationCancelledError,
        }
        exception_type = exception_types.get(self.status)
        if exception_type is not None:
            raise exception_type(self.failure_message or self.status.name, self)
