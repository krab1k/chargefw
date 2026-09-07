"""Read a molecular file, calculate an explicit charge model, and print the charges.

Run with: python docs/recipes/calculate_file.py molecule.sdf --format sdf \
    --method qeq --parameter-set QEq_original
"""

from __future__ import annotations

import argparse
from pathlib import Path

import chargefw


def calculate_file(
    input_path: str | Path,
    *,
    input_format: str,
    method: str,
    parameter_set: str,
) -> chargefw.CalculationResult:
    """Calculate full charges for every molecule in one supported input file."""

    molecules = chargefw.io.read(input_path, format=input_format)
    return chargefw.calculate(
        molecules,
        method=method,
        parameter_set=parameter_set,
        execution="full",
    )


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("input", type=Path)
    parser.add_argument("--format", required=True, choices=chargefw.io.INPUT_FORMATS)
    parser.add_argument("--method", required=True)
    parser.add_argument("--parameter-set", required=True)
    arguments = parser.parse_args()

    result = calculate_file(
        arguments.input,
        input_format=arguments.format,
        method=arguments.method,
        parameter_set=arguments.parameter_set,
    )
    for assignment in result.assignments:
        print(assignment.values)


if __name__ == "__main__":
    main()
