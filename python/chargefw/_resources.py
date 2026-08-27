"""Private package-resource discovery."""

import importlib.resources
from pathlib import Path


def default_parameter_directory() -> str:
    package_resources = importlib.resources.files("chargefw").joinpath("_data", "parameters")
    if package_resources.is_dir():
        return str(package_resources)
    return str(Path(__file__).resolve().parents[2] / "data" / "parameters")
