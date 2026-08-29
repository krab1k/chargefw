"""ChargeFW Python API."""

from ._chargefw import version as _native_version
from ._methods import (
    ExecutionIssue,
    ExecutionSupport,
    Method,
    MethodCatalog,
    MethodOption,
    MethodOptionCatalog,
    PrerequisiteIssue,
)
from ._parameters import ParameterSet, ParameterSetCatalog
from .calculation import (
    ApplicabilityReport,
    ApplicableMethod,
    Assessment,
    AssessmentReport,
    CalculationCancelledError,
    CalculationResult,
    CalculationTimings,
    ChargeFWError,
    EffectiveCalculation,
    ExecutionPolicy,
    InvalidInputError,
    NoExecutablePlanError,
    NumericalFailureError,
    RejectedMethod,
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
    "ExecutionSupport",
    "ExecutionPolicy",
    "ApplicableMethod",
    "RejectedMethod",
    "ApplicabilityReport",
    "EffectiveCalculation",
    "CalculationTimings",
    "AssessmentReport",
    "MethodOptionCatalog",
    "MethodOption",
    "MethodCatalog",
    "Method",
    "ParameterSetCatalog",
    "ParameterSet",
]
