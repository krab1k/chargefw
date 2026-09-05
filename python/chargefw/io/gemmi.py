"""Conversion from upstream Gemmi objects to ChargeFW molecules."""

from __future__ import annotations

from importlib import import_module
from typing import TYPE_CHECKING, Any

from ..core import MoleculeCollection
from . import BondStrategy, ConformerSelection, RecordSelection, parse_mmcif

if TYPE_CHECKING:
    import gemmi as _gemmi


def _require_gemmi() -> Any:
    try:
        return import_module("gemmi")
    except ModuleNotFoundError as error:
        if error.name != "gemmi":
            raise
        raise ImportError(
            "chargefw.io.gemmi requires the optional Gemmi Python integration; "
            "install it with `pip install chargefw[gemmi]`"
        ) from error


def from_structure(
    structure: _gemmi.Structure,
    *,
    source_name: str = "",
    selection: RecordSelection = "all",
    bonds: BondStrategy = "none",
    conformers: ConformerSelection = "all",
) -> MoleculeCollection:
    """Convert a ``gemmi.Structure`` to a one-molecule collection."""

    gemmi = _require_gemmi()
    if not isinstance(structure, gemmi.Structure):
        raise TypeError("structure must be a gemmi.Structure")
    return parse_mmcif(
        structure.make_mmcif_document().as_string(),
        source_name=source_name,
        selection=selection,
        bonds=bonds,
        conformers=conformers,
    )


def from_document(
    document: _gemmi.cif.Document,
    *,
    source_name: str = "",
    selection: RecordSelection = "all",
    bonds: BondStrategy = "none",
    conformers: ConformerSelection = "all",
) -> MoleculeCollection:
    """Convert a ``gemmi.cif.Document`` through ChargeFW's mmCIF parser."""

    gemmi = _require_gemmi()
    if not isinstance(document, gemmi.cif.Document):
        raise TypeError("document must be a gemmi.cif.Document")
    return parse_mmcif(
        document.as_string(),
        source_name=source_name,
        selection=selection,
        bonds=bonds,
        conformers=conformers,
    )


__all__ = ["from_structure", "from_document"]
