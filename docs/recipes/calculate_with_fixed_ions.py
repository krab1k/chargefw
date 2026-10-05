"""Calculate charges for a structure while treating selected ions as fixed sources.

Run with: python docs/recipes/calculate_with_fixed_ions.py structure.cif --format mmcif --ion MG
"""

import argparse
from pathlib import Path

import chargefw


def calculate_structure(
    input_path: str | Path,
    *,
    input_format: str,
    ions: list[str],
    method: str = "sqeqp",
    parameter_set: str | None = "SQEqp_Schindler2021_CCD_gen",
    selection: str = "all",
) -> chargefw.CalculationResult:
    """Read PDB/mmCIF input and calculate charges with the selected ions fixed."""

    molecules = chargefw.io.read(
        input_path,
        format=input_format,
        selection=selection,
        bonds="hybrid",
    )
    return chargefw.calculate(
        molecules,
        method=method,
        parameter_set=parameter_set,
        fixed_ions=ions,
    )


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("input", type=Path)
    parser.add_argument("--format", required=True, choices=("pdb", "mmcif"))
    parser.add_argument("--ion", action="append", required=True, dest="ions", metavar="ID")
    parser.add_argument("--method", default="sqeqp")
    parser.add_argument("--parameter-set", default="SQEqp_Schindler2021_CCD_gen")
    parser.add_argument(
        "--selection", choices=("all", "polymers-and-ligands", "polymers"), default="all"
    )
    parser.add_argument("--result-json", type=Path)
    arguments = parser.parse_args()

    result = calculate_structure(
        arguments.input,
        input_format=arguments.format,
        ions=arguments.ions,
        method=arguments.method,
        parameter_set=arguments.parameter_set,
        selection=arguments.selection,
    )
    if result.plan is None or result.plan.fixed_ions is None:
        raise RuntimeError("the calculation did not resolve any selected fixed ions")

    for warning in result.warnings:
        print(f"Warning: {warning.message}")
    for source in result.plan.fixed_ions.sources:
        print(
            f"Fixed ion: molecule={source.molecule_index} atom={source.atom_index} "
            f"charge={source.charge:+.1f} e"
        )
    for assignment in result.assignments:
        print(
            f"molecule={assignment.molecule_index} conformer={assignment.conformer_index}: "
            f"{assignment.values}"
        )
    if arguments.result_json is not None:
        chargefw.io.write(arguments.result_json, result, format="result-json")


if __name__ == "__main__":
    main()
