"""ChargeFW Python API."""

from ._chargefw import version as _native_version
from ._methods import (
    ExecutionIssue,
    Method,
    MethodCatalog,
    MethodOption,
    MethodOptionCatalog,
    PrerequisiteIssue,
)
from ._parameters import ParameterSet, ParameterSetCatalog
from .calculation import (
    Assessment,
    CalculationCancelledError,
    CalculationObserver,
    CalculationProgress,
    CalculationResult,
    CalculationTimings,
    ChargeFWError,
    ExecutedPlan,
    ExecutionPolicy,
    InvalidInputError,
    NoExecutablePlanError,
    NumericalFailureError,
    Plan,
    Rejection,
    RequestedCalculation,
    assess,
    calculate,
    methods,
    parameter_sets,
)
from .charges import ChargeAssignment
from .core import Molecule, MoleculeCollection, SourceIdentity

__version__ = _native_version()

__all__ = [
    "__version__",
    "Molecule",
    "MoleculeCollection",
    "SourceIdentity",
    "ChargeAssignment",
    "CalculationResult",
    "Assessment",
    "CalculationObserver",
    "CalculationProgress",
    "assess",
    "calculate",
    "methods",
    "parameter_sets",
    "RequestedCalculation",
    "ChargeFWError",
    "InvalidInputError",
    "NoExecutablePlanError",
    "NumericalFailureError",
    "CalculationCancelledError",
    "PrerequisiteIssue",
    "ExecutionIssue",
    "ExecutionPolicy",
    "Plan",
    "Rejection",
    "ExecutedPlan",
    "CalculationTimings",
    "MethodOptionCatalog",
    "MethodOption",
    "MethodCatalog",
    "Method",
    "ParameterSetCatalog",
    "ParameterSet",
]
