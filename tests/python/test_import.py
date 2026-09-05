"""Focused smoke test for the package skeleton."""

import os
from pathlib import Path

import chargefw
import chargefw.calculation
import chargefw.core
import chargefw.io
import chargefw.io.gemmi
from chargefw import _chargefw

expected_version = os.environ.get("CHARGEFW_EXPECTED_VERSION")


def test_import_surface() -> None:
    assert isinstance(chargefw.__version__, str)
    assert chargefw.__version__
    if expected_version is not None:
        assert chargefw.__version__ == expected_version
    assert chargefw.__all__ == [
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
        "io",
    ]
    assert chargefw.Molecule is chargefw.core.Molecule
    assert chargefw.RequestedCalculation is chargefw.calculation.RequestedCalculation
    assert chargefw.CalculationObserver is chargefw.calculation.CalculationObserver
    assert chargefw.CalculationProgress is chargefw.calculation.CalculationProgress
    assert isinstance(chargefw.methods, chargefw.MethodCatalog)
    assert isinstance(chargefw.parameter_sets, chargefw.ParameterSetCatalog)
    assert chargefw.calculate is chargefw.calculation.calculate
    assert chargefw.assess is chargefw.calculation.assess
    assert chargefw.io is chargefw.io
    assert chargefw.io.__all__ == [
        "InputFormat",
        "RecordSelection",
        "BondStrategy",
        "ConformerSelection",
        "parse",
        "read",
    ]
    assert chargefw.io.gemmi.__all__ == ["from_structure", "from_document"]
    assert not hasattr(_chargefw.calculation, "ExecutionMode")
    assert not hasattr(_chargefw.methods, "MethodOptionType")
    assert not hasattr(_chargefw.adapters, "BondStrategy")
    assert (Path(chargefw.__file__).parent / "_chargefw" / "__init__.pyi").is_file()
    assert (Path(chargefw.__file__).parent / "_chargefw" / "adapters.pyi").is_file()
    assert (Path(chargefw.__file__).parent / "py.typed").is_file()


if __name__ == "__main__":
    test_import_surface()
