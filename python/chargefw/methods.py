"""Method applicability and execution report values."""

from __future__ import annotations

from dataclasses import dataclass

from ._chargefw import calculation as _native_calculation
from ._chargefw import methods as _native_methods

PrerequisiteIssueKind = _native_methods.PrerequisiteIssueKind
ExecutionAvailability = _native_methods.ExecutionAvailability
ExecutionIssueKind = _native_methods.ExecutionIssueKind
ExecutionMode = _native_calculation.ExecutionMode
for _enum in (PrerequisiteIssueKind, ExecutionAvailability, ExecutionIssueKind):
    _enum.__module__ = __name__


@dataclass(frozen=True, slots=True)
class PrerequisiteIssue:
    kind: PrerequisiteIssueKind
    message: str
    molecule_index: int | None = None
    atom_index: int | None = None
    bond_index: int | None = None
    conformer_index: int | None = None


@dataclass(frozen=True, slots=True)
class ExecutionIssue:
    kind: ExecutionIssueKind
    message: str
    molecule_index: int | None = None


@dataclass(frozen=True, slots=True)
class ExecutionAssessment:
    mode: ExecutionMode
    availability: ExecutionAvailability
    issues: tuple[ExecutionIssue, ...] = ()
