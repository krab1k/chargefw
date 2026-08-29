"""Focused smoke test for the package skeleton."""

import os
from pathlib import Path

import chargefw
import chargefw.adapters
import chargefw.adapters.gemmi
import chargefw.calculation
import chargefw.core

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
        "CalculationOptions",
        "ChargeAssignment",
        "CalculationResult",
        "Assessment",
        "assess",
        "calculate",
        "methods",
        "parameter_sets",
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
        "MethodOptionType",
        "MethodOptionCatalog",
        "MethodOptionDescriptor",
        "MethodCatalog",
        "MethodDescriptor",
        "ParameterSetCatalog",
        "ParameterSetDescriptor",
    ]
    assert chargefw.Molecule is chargefw.core.Molecule
    assert chargefw.CalculationOptions is chargefw.calculation.CalculationOptions
    assert isinstance(chargefw.methods, chargefw.MethodCatalog)
    assert isinstance(chargefw.parameter_sets, chargefw.ParameterSetCatalog)
    assert chargefw.calculate is chargefw.calculation.calculate
    assert chargefw.assess is chargefw.calculation.assess
    assert chargefw.ExecutionMode.__module__ == "chargefw.calculation"
    assert chargefw.PrerequisiteIssueKind.__module__ == "chargefw"
    assert chargefw.MethodOptionType.__module__ == "chargefw"
    assert chargefw.adapters.read_pdb is chargefw.adapters.gemmi.read_pdb
    assert chargefw.adapters.BondStrategy.__module__ == "chargefw.adapters.gemmi"
    assert chargefw.adapters.__all__ == [
        "RecordSelection",
        "BondStrategy",
        "ConformerSelection",
        "read_pdb_string",
        "read_pdb",
        "read_mmcif_string",
        "read_mmcif",
        "from_structure",
        "from_document",
    ]
    assert (Path(chargefw.__file__).parent / "_chargefw" / "__init__.pyi").is_file()
    assert (Path(chargefw.__file__).parent / "_chargefw" / "adapters.pyi").is_file()
    assert (Path(chargefw.__file__).parent / "py.typed").is_file()


if __name__ == "__main__":
    test_import_surface()
