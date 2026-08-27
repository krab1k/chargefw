"""ChargeFW Python API."""

from ._chargefw import version as _native_version
from .calculation import (
    ApplicableCandidate,
    ApplicabilityReport,
    Assessment,
    AssessmentReport,
    CalculationCancelledError,
    CalculationOptions,
    CalculationResult,
    CalculationTimings,
    Calculator,
    ChargeCorrectionPolicy,
    ChargeFWError,
    EffectiveCalculation,
    ExecutionMode,
    ExecutionPolicy,
    ExecutionSelectionKind,
    ExecutionStatus,
    InvalidInputError,
    NoExecutablePlanError,
    NumericalFailureError,
    RejectedCandidate,
)
from .charges import ChargeAssignment
from .core import Molecule, MoleculeCollection, SourceIdentity
from .methods import (
    ExecutionAssessment,
    ExecutionAvailability,
    ExecutionIssue,
    ExecutionIssueKind,
    PrerequisiteIssue,
    PrerequisiteIssueKind,
)

__version__ = _native_version()

__all__ = [
    "__version__",
    "Molecule",
    "MoleculeCollection",
    "SourceIdentity",
    "CalculationOptions",
    "ChargeAssignment",
    "CalculationResult",
    "Assessment",
    "Calculator",
    "ChargeFWError",
    "InvalidInputError",
    "NoExecutablePlanError",
    "NumericalFailureError",
    "CalculationCancelledError",
    "ExecutionSelectionKind",
    "ExecutionMode",
    "ChargeCorrectionPolicy",
    "ExecutionStatus",
    "PrerequisiteIssueKind",
    "ExecutionAvailability",
    "ExecutionIssueKind",
    "PrerequisiteIssue",
    "ExecutionIssue",
    "ExecutionAssessment",
    "ExecutionPolicy",
    "ApplicableCandidate",
    "RejectedCandidate",
    "ApplicabilityReport",
    "EffectiveCalculation",
    "CalculationTimings",
    "AssessmentReport",
]
