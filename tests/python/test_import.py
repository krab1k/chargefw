"""Focused smoke test for the package skeleton."""

import os
from pathlib import Path

import chargefw
import chargefw.calculation
import chargefw.core
import chargefw.io
import chargefw.io.gemmi
import chargefw.io.rdkit

expected_version = os.environ.get("CHARGEFW_EXPECTED_VERSION")


def test_import_surface() -> None:
    assert isinstance(chargefw.__version__, str)
    assert chargefw.__version__
    if expected_version is not None:
        assert chargefw.__version__ == expected_version
    assert chargefw.Molecule is chargefw.core.Molecule
    assert chargefw.RequestedCalculation is chargefw.calculation.RequestedCalculation
    assert chargefw.CalculationObserver is chargefw.calculation.CalculationObserver
    assert chargefw.CalculationProgress is chargefw.calculation.CalculationProgress
    assert chargefw.calculate is chargefw.calculation.calculate
    assert chargefw.assess is chargefw.calculation.assess
    assert (Path(chargefw.__file__).parent / "_chargefw" / "__init__.pyi").is_file()
    assert (Path(chargefw.__file__).parent / "_chargefw" / "adapters.pyi").is_file()
    assert (Path(chargefw.__file__).parent / "py.typed").is_file()


if __name__ == "__main__":
    test_import_surface()
