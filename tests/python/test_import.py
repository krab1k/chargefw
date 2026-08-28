"""Focused smoke test for the package skeleton."""

import os
from pathlib import Path

import chargefw
import chargefw.adapters
import chargefw.adapters.gemmi
import chargefw.calculation
import chargefw.core
import chargefw.methods
import chargefw.parameters

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
        "MethodOptionType",
        "MethodOptionDescriptor",
        "MethodDescriptor",
        "method_descriptors",
        "ParameterSet",
        "ParameterSetDescriptor",
        "load_parameter_set",
        "load_parameter_sets",
    ]
    assert chargefw.Molecule is chargefw.core.Molecule
    assert chargefw.CalculationOptions is chargefw.calculation.CalculationOptions
    assert chargefw.ExecutionIssue is chargefw.methods.ExecutionIssue
    assert chargefw.ExecutionMode.__module__ == "chargefw.calculation"
    assert chargefw.PrerequisiteIssueKind.__module__ == "chargefw.methods"
    assert chargefw.MethodOptionType.__module__ == "chargefw.methods"
    assert chargefw.ParameterSet is chargefw.parameters.ParameterSet
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
