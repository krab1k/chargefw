"""Typed private payloads returned by the nanobind extension."""

from __future__ import annotations

from typing import TypedDict

import numpy as np
from numpy.typing import NDArray

from ._types import (
    ChargeCorrection,
    ExecutionIssueKind,
    ExecutionMode,
    ExecutionStatus,
    MethodOptionType,
    MethodOptionValue,
    PrerequisiteIssueKind,
)


class PrerequisiteIssuePayload(TypedDict):
    kind: PrerequisiteIssueKind
    message: str
    molecule_index: int | None
    atom_index: int | None
    bond_index: int | None
    conformer_index: int | None


class ExecutionIssuePayload(TypedDict):
    kind: ExecutionIssueKind
    message: str
    molecule_index: int | None


class ExecutionPolicyPayload(TypedDict):
    mode: ExecutionMode
    radius: float | None
    charge_correction: ChargeCorrection


class EffectiveCalculationPayload(TypedDict):
    method_id: str
    parameter_set_id: str | None
    method_options: dict[str, MethodOptionValue]
    execution_policy: ExecutionPolicyPayload
    execution_issues: list[ExecutionIssuePayload]


class ExecutionPlanPayload(TypedDict):
    method_id: str
    parameter_set_id: str | None
    method_options: dict[str, MethodOptionValue]
    execution_policy: ExecutionPolicyPayload
    warnings: list[ExecutionIssuePayload]


class RejectionPayload(TypedDict):
    method_id: str
    parameter_set_id: str | None
    execution_policy: ExecutionPolicyPayload | None
    prerequisite_issues: list[PrerequisiteIssuePayload]
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
    status: ExecutionStatus
    rejections: list[RejectionPayload]
    failure_message: str | None
    metrics: CalculationMetricsPayload
    effective: EffectiveCalculationPayload | None
    charges: ChargeSetPayload | None


class AssessmentReportPayload(TypedDict):
    rejections: list[RejectionPayload]
    applicability_seconds: float


class MethodOptionDescriptorPayload(TypedDict):
    id: str
    description: str
    type: MethodOptionType
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
