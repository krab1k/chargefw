"""Focused smoke test for the package skeleton."""

import sys

import chargefw


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
]
