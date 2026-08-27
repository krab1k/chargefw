"""Focused smoke test for the package skeleton."""

import sys

import chargefw
import chargefw.calculation
import chargefw.core
import chargefw.methods
import chargefw.parameters


expected_version = sys.argv[1]

assert chargefw.__version__ == expected_version
assert chargefw.__all__ == [
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
assert chargefw.Molecule is chargefw.core.Molecule
assert chargefw.CalculationOptions is chargefw.calculation.CalculationOptions
assert chargefw.ExecutionIssue is chargefw.methods.ExecutionIssue
assert chargefw.ExecutionMode.__module__ == "chargefw.calculation"
assert chargefw.PrerequisiteIssueKind.__module__ == "chargefw.methods"
