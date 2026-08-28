"""Typed private payloads returned by the nanobind extension."""

from __future__ import annotations

from typing import TypedDict

import numpy as np
from numpy.typing import NDArray

from ._chargefw import calculation as _native_calculation
from ._chargefw import methods as _native_methods
from ._types import MethodOptionValue


class PrerequisiteIssuePayload(TypedDict):
    kind: _native_methods.PrerequisiteIssueKind
    message: str
    molecule_index: int | None
    atom_index: int | None
    bond_index: int | None
    conformer_index: int | None


class ExecutionIssuePayload(TypedDict):
    kind: _native_methods.ExecutionIssueKind
    message: str
    molecule_index: int | None


class ExecutionAssessmentPayload(TypedDict):
    mode: _native_calculation.ExecutionMode
    availability: _native_methods.ExecutionAvailability
    issues: list[ExecutionIssuePayload]


class ApplicableCandidatePayload(TypedDict):
    method_id: str
    parameter_set_id: str | None
    execution_assessments: list[ExecutionAssessmentPayload]


class RejectedCandidatePayload(TypedDict):
    method_id: str
    parameter_set_id: str | None
    issues: list[PrerequisiteIssuePayload]


class ApplicabilityReportPayload(TypedDict):
    applicable: list[ApplicableCandidatePayload]
    rejected: list[RejectedCandidatePayload]
    selected_candidate_index: int | None


class ExecutionPolicyPayload(TypedDict):
    mode: _native_calculation.ExecutionMode
    radius: float | None
    charge_correction: _native_calculation.ChargeCorrectionPolicy


class EffectiveCalculationPayload(TypedDict):
    method_id: str
    parameter_set_id: str | None
    method_options: dict[str, MethodOptionValue]
    execution_policy: ExecutionPolicyPayload
    execution_issues: list[ExecutionIssuePayload]


class ChargeAssignmentPayload(TypedDict):
    molecule_index: int
    conformer_index: int | None
    values: NDArray[np.float64]


class ChargeSetPayload(TypedDict):
    method_id: str
    parameter_set_id: str | None
    assignments: list[ChargeAssignmentPayload]


class CalculationMetricsPayload(TypedDict):
    applicability_seconds: float
    computation_seconds: float


class ExecutionResultPayload(TypedDict):
    status: _native_calculation.ExecutionStatus
    applicability: ApplicabilityReportPayload
    failure_message: str | None
    metrics: CalculationMetricsPayload
    effective: EffectiveCalculationPayload | None
    charges: ChargeSetPayload | None


class AssessmentReportPayload(TypedDict):
    applicability: ApplicabilityReportPayload
    execution_policy: ExecutionPolicyPayload | None
    execution_issues: list[ExecutionIssuePayload]
    applicability_seconds: float
    executable: bool


class MethodOptionDescriptorPayload(TypedDict):
    id: str
    description: str
    type: _native_methods.MethodOptionType
    default: MethodOptionValue
    choices: list[MethodOptionValue]
    minimum: MethodOptionValue | None
    minimum_inclusive: bool
    maximum: MethodOptionValue | None
    maximum_inclusive: bool


class MethodDescriptorPayload(TypedDict):
    id: str
    name: str
    full_name: str
    publication: str | None
    priority: int
    requires_coordinates: bool
    supports_cutoff: bool
    supports_cover: bool
    options: list[MethodOptionDescriptorPayload]


class ParameterSetDescriptorPayload(TypedDict):
    id: str
    method_id: str
    name: str
    publication: str
    notes: str
    priority: int
