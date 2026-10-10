"""Focused smoke test for the installed package."""

import os
from pathlib import Path

import chargefw
import chargefw.io.gemmi  # noqa: F401
import chargefw.io.rdkit  # noqa: F401


def test_version_and_typing_marker() -> None:
    assert isinstance(chargefw.__version__, str)
    assert chargefw.__version__
    expected_version = os.environ.get("CHARGEFW_EXPECTED_VERSION")
    if expected_version is not None:
        assert chargefw.__version__ == expected_version
    assert (Path(chargefw.__file__).parent / "py.typed").is_file()
