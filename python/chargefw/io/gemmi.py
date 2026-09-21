"""Conversion from upstream Gemmi objects to ChargeFW molecules."""

from __future__ import annotations

from importlib import import_module
from typing import TYPE_CHECKING, Any, cast

from .._chargefw import adapters as _native_adapters
from ..core import MoleculeCollection
from . import (
    BondStrategy,
    ConformerSelection,
    RecordSelection,
    dumps,
    parse,
)

if TYPE_CHECKING:
    import gemmi as _gemmi

    from ..calculation import CalculationResult


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
    return parse(
        structure.make_mmcif_document().as_string(),
        format="mmcif",
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
    return parse(
        document.as_string(),
        format="mmcif",
        source_name=source_name,
        selection=selection,
        bonds=bonds,
        conformers=conformers,
    )


def to_document(result: CalculationResult) -> _gemmi.cif.Document:
    """Create a fresh Gemmi mmCIF document from a calculation result."""

    from ..calculation import CalculationResult

    gemmi = _require_gemmi()
    if not isinstance(result, CalculationResult):
        raise TypeError("result must be a CalculationResult")
    return cast("_gemmi.cif.Document", gemmi.cif.read_string(dumps(result, format="mmcif")))


def attach_charges(
    document: _gemmi.cif.Document,
    result: CalculationResult,
    *,
    overwrite: bool = False,
) -> None:
    """Attach charges to the unchanged Gemmi document used for the calculation input."""

    from ..calculation import CalculationResult

    gemmi = _require_gemmi()
    if not isinstance(document, gemmi.cif.Document):
        raise TypeError("document must be a gemmi.cif.Document")
    if not isinstance(result, CalculationResult):
        raise TypeError("result must be a CalculationResult")
    if not isinstance(overwrite, bool):
        raise TypeError("overwrite must be a bool")
    attached = gemmi.cif.read_string(
        _native_adapters._attach_mmcif(document.as_string(), result._native, overwrite)
    )
    document.clear()
    for block in attached:
        document.add_copied_block(block)


__all__ = ["from_structure", "from_document", "to_document", "attach_charges"]
