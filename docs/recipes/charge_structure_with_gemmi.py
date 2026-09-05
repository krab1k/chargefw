"""Calculate charges for a Gemmi structure and write a charged mmCIF document.

Run with: python docs/recipes/charge_structure_with_gemmi.py input.cif charged.cif
"""

from __future__ import annotations

import argparse
from pathlib import Path

import chargefw
import gemmi
from chargefw.io import gemmi as chargefw_gemmi


def charge_structure(
    input_path: str | Path, output_path: str | Path, *, method: str = "eem"
) -> None:
    """Load a structure, calculate charges, and write an enriched mmCIF document."""

    input_path = Path(input_path)
    structure = gemmi.read_structure(str(input_path))
    document = structure.make_mmcif_document()
    molecules = chargefw_gemmi.from_document(document, source_name=str(input_path), bonds="hybrid")
    result = chargefw.calculate(molecules, method=method, execution="full")
    chargefw_gemmi.attach_charges(document, result)
    document.write_file(str(output_path))


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("input", type=Path)
    parser.add_argument("output", type=Path)
    parser.add_argument("--method", default="eem")
    arguments = parser.parse_args()
    charge_structure(arguments.input, arguments.output, method=arguments.method)


if __name__ == "__main__":
    main()
