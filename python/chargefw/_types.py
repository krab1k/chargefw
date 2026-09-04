"""Shared private type aliases."""

from typing import Literal, TypeAlias

MethodOptionValue: TypeAlias = bool | int | float | str
Execution: TypeAlias = Literal["auto", "full", "cutoff", "cover"]
ChargeCorrection: TypeAlias = Literal["none", "uniform"]
ParameterMatching: TypeAlias = Literal["strict", "permissive"]
ExecutionMode: TypeAlias = Literal["full", "cutoff", "cover"]
CalculationPhase: TypeAlias = Literal[
    "computation_started",
    "computation_finished",
    "target_started",
    "target_finished",
    "fragment_progress",
]
ExecutionIssueKind: TypeAlias = Literal["resource_threshold_exceeded", "unsupported_execution_mode"]
PrerequisiteIssueKind: TypeAlias = Literal[
    "invalid_options",
    "missing_feature",
    "invalid_geometry",
    "unsupported_molecule",
    "missing_parameters",
    "parameter_classification_failed",
]
MethodOptionType: TypeAlias = Literal["boolean", "integer", "floating_point", "string"]
ExecutionStatus: TypeAlias = Literal[
    "success",
    "invalid_input_or_request",
    "no_executable_plan",
    "numerical_failure",
    "cancelled",
]
