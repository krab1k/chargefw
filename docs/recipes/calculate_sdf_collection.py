"""Calculate EEM charges for every record in a small SDF collection.

Run with: python docs/recipes/calculate_sdf_collection.py input.sdf result.json
"""

import argparse
from pathlib import Path

import chargefw


def calculate_sdf_collection(input_path: str | Path, output_path: str | Path) -> None:
    """Calculate EEM charges for the collection and write complete result JSON.

    One parameter set and execution mode are selected automatically for all records.
    """

    molecules = chargefw.io.read(input_path, format="sdf")
    result = chargefw.calculate(
        molecules,
        method="eem",
    )
    chargefw.io.write(output_path, result, format="result-json")
    print(
        f"Calculated {len(result.assignments)} assignment(s) for {len(molecules)} "
        f"molecule(s) with {result.plan.method.id}"
    )


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("input", type=Path)
    parser.add_argument("output", type=Path)
    arguments = parser.parse_args()
    calculate_sdf_collection(arguments.input, arguments.output)


if __name__ == "__main__":
    main()
