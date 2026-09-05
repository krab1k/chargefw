"""Native molecular input and generated calculation output."""

from __future__ import annotations

from operator import index as as_index
from os import PathLike
from pathlib import Path
from typing import TYPE_CHECKING, Literal, TypeAlias

from .._chargefw import adapters as _native_adapters
from ..core import Molecule, MoleculeCollection

if TYPE_CHECKING:
    from ..calculation import CalculationResult

InputFormat: TypeAlias = Literal["mol", "sdf", "mol2", "molecule-json", "pdb", "mmcif"]
OutputFormat: TypeAlias = Literal["sdf", "mol2", "mmcif", "result-json"]
SdfVersion: TypeAlias = Literal["v2000", "v3000"]
RecordSelection: TypeAlias = Literal["all", "polymers-and-ligands", "polymers"]
BondStrategy: TypeAlias = Literal["none", "explicit", "templates", "hybrid"]
ConformerSelection: TypeAlias = Literal["first", "all"]

_INPUT_FORMATS = frozenset(("mol", "sdf", "mol2", "molecule-json", "pdb", "mmcif"))
_STRUCTURAL_FORMATS = frozenset(("pdb", "mmcif"))
_MULTI_CONFORMER_FORMATS = frozenset(("molecule-json", "pdb", "mmcif"))
_OUTPUT_FORMATS = frozenset(("sdf", "mol2", "mmcif", "result-json"))


def _molecule(payload: _native_adapters.MoleculePayload) -> Molecule:
    return Molecule(
        atomic_numbers=payload["atomic_numbers"],
        formal_charges=payload["formal_charges"],
        bonds=payload["bonds"],
        coordinates=payload["coordinates"],
        name=payload["name"],
        atom_names=payload["atom_names"],
        conformer_names=payload["conformer_names"],
        source_name=payload["source"],
        record_index=payload["record_index"],
        record_id=payload["record_id"],
    )


def _collection(
    payloads: list[_native_adapters.MoleculePayload], source_name: str
) -> MoleculeCollection:
    return MoleculeCollection((_molecule(payload) for payload in payloads), name=source_name)


def _validate_options(
    format: InputFormat,
    selection: RecordSelection,
    bonds: BondStrategy,
    conformers: ConformerSelection,
) -> None:
    if not isinstance(format, str):
        raise TypeError("format must be a string")
    if format not in _INPUT_FORMATS:
        raise ValueError(f"unsupported molecular input format: {format}")
    if format not in _STRUCTURAL_FORMATS:
        if selection != "all":
            raise ValueError("selection is only supported for PDB and mmCIF input")
        if bonds != "none":
            raise ValueError("bonds is only supported for PDB and mmCIF input")
    if format not in _MULTI_CONFORMER_FORMATS and conformers != "all":
        raise ValueError("conformers is only supported for molecule JSON, PDB, and mmCIF input")


def parse(
    contents: str,
    *,
    format: InputFormat,
    source_name: str = "",
    selection: RecordSelection = "all",
    bonds: BondStrategy = "none",
    conformers: ConformerSelection = "all",
) -> MoleculeCollection:
    """Parse molecular text using an explicitly selected native reader."""

    if not isinstance(contents, str):
        raise TypeError("contents must be a string")
    if not isinstance(source_name, str):
        raise TypeError("source_name must be a string")
    _validate_options(format, selection, bonds, conformers)
    return _collection(
        _native_adapters._parse(contents, source_name, format, selection, bonds, conformers),
        source_name,
    )


def read(
    path: str | PathLike[str],
    *,
    format: InputFormat,
    selection: RecordSelection = "all",
    bonds: BondStrategy = "none",
    conformers: ConformerSelection = "all",
) -> MoleculeCollection:
    """Read UTF-8 molecular text using an explicitly selected native reader."""

    if not isinstance(path, (str, PathLike)):
        raise TypeError("path must be a string or path-like value")
    source_path = Path(path)
    return parse(
        source_path.read_text(encoding="utf-8"),
        format=format,
        source_name=str(source_path),
        selection=selection,
        bonds=bonds,
        conformers=conformers,
    )


def dumps(
    result: CalculationResult,
    *,
    format: OutputFormat,
    conformer: int | None = None,
    sdf_version: SdfVersion | None = None,
) -> str:
    """Serialize a calculation result through a native generated-output writer."""

    from ..calculation import CalculationResult

    if not isinstance(result, CalculationResult):
        raise TypeError("result must be a CalculationResult")
    if not isinstance(format, str):
        raise TypeError("format must be a string")
    if format not in _OUTPUT_FORMATS:
        raise ValueError(f"unsupported calculation output format: {format}")
    if isinstance(conformer, bool):
        raise TypeError("conformer must be an integer or None")
    if conformer is not None:
        try:
            conformer = as_index(conformer)
        except TypeError as error:
            raise TypeError("conformer must be an integer or None") from error
        if conformer < 0:
            raise ValueError("conformer must be non-negative")
    if conformer is not None and format not in ("sdf", "mol2"):
        raise ValueError("conformer is only supported for SDF and MOL2 output")
    if sdf_version is not None and format != "sdf":
        raise ValueError("sdf_version is only supported for SDF output")
    if sdf_version is not None and sdf_version not in ("v2000", "v3000"):
        raise ValueError("sdf_version must be 'v2000', 'v3000', or None")
    identities: list[tuple[str, int, str]] = []
    for molecule in result.molecules:
        if molecule.record_id is not None and not isinstance(molecule.record_id, str):
            raise TypeError("serialized molecule record IDs must be strings or None")
        identities.append((molecule.source_name, molecule.record_index, molecule.record_id or ""))
    return _native_adapters._dumps(
        result._native,
        result.molecules._native_molecules,
        identities,
        result._requested_payload,
        format,
        conformer,
        sdf_version or "v3000",
    )


def write(
    path: str | PathLike[str],
    result: CalculationResult,
    *,
    format: OutputFormat,
    conformer: int | None = None,
    sdf_version: SdfVersion | None = None,
) -> None:
    """Serialize a calculation result to a UTF-8 text file."""

    if not isinstance(path, (str, PathLike)):
        raise TypeError("path must be a string or path-like value")
    Path(path).write_text(
        dumps(
            result,
            format=format,
            conformer=conformer,
            sdf_version=sdf_version,
        ),
        encoding="utf-8",
    )


__all__ = [
    "InputFormat",
    "OutputFormat",
    "SdfVersion",
    "RecordSelection",
    "BondStrategy",
    "ConformerSelection",
    "parse",
    "read",
    "dumps",
    "write",
]
