"""Conversion from upstream Gemmi objects to ChargeFW molecules."""

from __future__ import annotations

import gemmi as _gemmi

from ..core import MoleculeCollection
from . import BondStrategy, ConformerSelection, RecordSelection, parse_mmcif


def from_structure(
    structure: _gemmi.Structure,
    *,
    source_name: str = "",
    selection: RecordSelection = "all",
    bonds: BondStrategy = "none",
    conformers: ConformerSelection = "all",
) -> MoleculeCollection:
    """Convert a ``gemmi.Structure`` to a one-molecule collection."""

    if not isinstance(structure, _gemmi.Structure):
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

    if not isinstance(document, _gemmi.cif.Document):
        raise TypeError("document must be a gemmi.cif.Document")
    return parse_mmcif(
        document.as_string(),
        source_name=source_name,
        selection=selection,
        bonds=bonds,
        conformers=conformers,
    )


__all__ = ["from_structure", "from_document"]
