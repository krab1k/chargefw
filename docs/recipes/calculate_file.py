"""Read a molecular file, calculate charges, and print them per atom.

Run with: python docs/recipes/calculate_file.py molecules.sdf --format sdf --method eem
"""

import argparse
from pathlib import Path

import chargefw


def calculate_file(
    input_path: str | Path,
    *,
    input_format: str,
    method: str,
    parameter_set: str | None = None,
) -> tuple[chargefw.MoleculeCollection, chargefw.CalculationResult]:
    """Calculate charges for every molecule in one supported input file.

    One method, parameter set, and execution mode apply to the whole collection. The parameter set
    and execution mode are selected automatically when omitted.
    """

    molecules = chargefw.io.read(input_path, format=input_format)
    result = chargefw.calculate(
        molecules,
        method=method,
        parameter_set=parameter_set,
    )
    return molecules, result


def print_charges(
    molecules: chargefw.MoleculeCollection, result: chargefw.CalculationResult
) -> None:
    """Print one row per atom and assignment in source order."""

    print("molecule  conformer  atom_id  atom_name  element      charge")
    for molecule, assignments in zip(molecules, result.assignments_by_molecule):
        for assignment in assignments:
            # Geometry-independent methods assign charges to the molecule, not to one conformer.
            conformer = "-" if assignment.conformer_index is None else assignment.conformer_index
            for atom_index, charge in enumerate(assignment.values):
                element = chargefw.periodic_table[molecule.atomic_numbers[atom_index]].symbol
                name = molecule.atom_names[atom_index] if molecule.atom_names else ""
                print(
                    f"{assignment.molecule_index:>8}  {conformer:>9}  "
                    f"{molecule.atom_ids[atom_index]!s:>7}  {name or '-':<9}  "
                    f"{element:<7}  {charge:>10.6f}"
                )


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("input", type=Path)
    parser.add_argument("--format", required=True, choices=chargefw.io.INPUT_FORMATS)
    parser.add_argument("--method", required=True)
    parser.add_argument("--parameter-set")
    parser.add_argument("--result-json", type=Path, help="Also write complete result JSON")
    arguments = parser.parse_args()

    molecules, result = calculate_file(
        arguments.input,
        input_format=arguments.format,
        method=arguments.method,
        parameter_set=arguments.parameter_set,
    )
    print_charges(molecules, result)
    if arguments.result_json is not None:
        chargefw.io.write(arguments.result_json, result, format="result-json")


if __name__ == "__main__":
    main()
