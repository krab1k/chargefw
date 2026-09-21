"""Native molecular input and generated calculation output."""

from __future__ import annotations

from dataclasses import dataclass
from os import PathLike
from pathlib import Path
from typing import TYPE_CHECKING, Literal, TypeAlias

from .._chargefw import adapters as _native_adapters
from ..core import (
    Molecule,
    MoleculeCollection,
    SourceAtomReference,
    SourceConformerReference,
    SourceMapping,
)

if TYPE_CHECKING:
    from ..calculation import CalculationResult

InputFormat: TypeAlias = Literal["mol", "sdf", "mol2", "molecule-json", "pdb", "mmcif"]
OutputFormat: TypeAlias = Literal["sdf", "mol2", "mmcif", "result-json"]
SdfVersion: TypeAlias = Literal["v2000", "v3000"]
RecordSelection: TypeAlias = Literal["all", "polymers-and-ligands", "polymers"]
BondStrategy: TypeAlias = Literal["none", "explicit", "templates", "hybrid"]
ConformerSelection: TypeAlias = Literal["first", "all"]

INPUT_FORMATS: tuple[InputFormat, ...] = ("mol", "sdf", "mol2", "molecule-json", "pdb", "mmcif")
OUTPUT_FORMATS: tuple[OutputFormat, ...] = ("sdf", "mol2", "mmcif", "result-json")
_STRUCTURAL_FORMATS = frozenset(("pdb", "mmcif"))
_MULTI_CONFORMER_FORMATS = frozenset(("molecule-json", "pdb", "mmcif"))


@dataclass(frozen=True, slots=True)
class _InputMetadata:
    diagnostics: tuple[tuple[str, str, int | None], ...]
    structural_input: tuple[RecordSelection, BondStrategy] | None
    conformers: ConformerSelection


class _ImportedMolecule(Molecule):
    __slots__ = ("_input_metadata",)

    _input_metadata: _InputMetadata


def _molecule(
    payload: _native_adapters.MoleculePayload,
    structural_input: tuple[RecordSelection, BondStrategy] | None,
    conformers: ConformerSelection,
) -> Molecule:
    coordinates = payload["coordinates"] or None
    import_metadata = payload["import_metadata"]
    source_mapping = (
        None
        if import_metadata is None
        else SourceMapping(
            format=import_metadata["format"],
            atoms=tuple(
                SourceAtomReference(position, source_id)
                for position, source_id in import_metadata["atoms"]
            ),
            conformers=tuple(
                SourceConformerReference(
                    position,
                    source_id,
                    tuple(
                        SourceAtomReference(site_position, site_id)
                        for site_position, site_id in sites
                    ),
                )
                for position, source_id, sites in import_metadata["conformers"]
            ),
            source_connectivity=import_metadata["source_connectivity"],
            record_selection=import_metadata["record_selection"],
            conformer_selection=import_metadata["conformer_selection"],
            bond_strategy=import_metadata["bond_strategy"],
        )
    )
    result = _ImportedMolecule(
        atomic_numbers=payload["atomic_numbers"],
        formal_charges=payload["formal_charges"],
        bonds=payload["bonds"],
        coordinates=coordinates,
        name=payload["name"],
        atom_names=payload["atom_names"],
        conformer_names=payload["conformer_names"],
        source_name=payload["source"],
        record_index=payload["record_index"],
        record_id=payload["record_id"],
        atom_ids=(
            None
            if source_mapping is None
            else tuple(
                reference.position if reference.id is None else reference.id
                for reference in source_mapping.atoms
            )
        ),
    )
    object.__setattr__(result, "_source_mapping", source_mapping)
    object.__setattr__(
        result,
        "_input_metadata",
        _InputMetadata(
            diagnostics=tuple(payload["diagnostics"]),
            structural_input=structural_input,
            conformers=conformers,
        ),
    )
    return result


def _collection(
    payloads: list[_native_adapters.MoleculePayload],
    source_name: str,
    structural_input: tuple[RecordSelection, BondStrategy] | None,
    conformers: ConformerSelection,
) -> MoleculeCollection:
    molecules = tuple(_molecule(payload, structural_input, conformers) for payload in payloads)
    return MoleculeCollection(molecules, source_name)


def _validate_options(
    format: InputFormat,
    selection: RecordSelection,
    bonds: BondStrategy,
    conformers: ConformerSelection,
) -> None:
    if not isinstance(format, str):
        raise TypeError("format must be a string")
    if format not in INPUT_FORMATS:
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
    payloads = _native_adapters._parse(contents, source_name, format, selection, bonds, conformers)
    structural_input = (selection, bonds) if format in _STRUCTURAL_FORMATS else None
    return _collection(payloads, source_name, structural_input, conformers)


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
    sdf_version: SdfVersion | None = None,
) -> str:
    """Serialize a calculation result through a native generated-output writer."""

    from ..calculation import CalculationResult

    if not isinstance(result, CalculationResult):
        raise TypeError("result must be a CalculationResult")
    if not isinstance(format, str):
        raise TypeError("format must be a string")
    if format not in OUTPUT_FORMATS:
        raise ValueError(f"unsupported calculation output format: {format}")
    if sdf_version is not None and format != "sdf":
        raise ValueError("sdf_version is only supported for SDF output")
    if sdf_version is not None and sdf_version not in ("v2000", "v3000"):
        raise ValueError("sdf_version must be 'v2000', 'v3000', or None")
    identities: list[tuple[str, int, str]] = []
    for molecule in result.molecules:
        record_id = molecule.record_id
        if format == "result-json" and record_id is not None and not isinstance(record_id, str):
            raise TypeError("result JSON record IDs must be strings or None")
        serialized_id = record_id if isinstance(record_id, str) else ""
        identities.append((molecule.source_name, molecule.record_index, serialized_id))
    metadata = tuple(
        molecule._input_metadata if isinstance(molecule, _ImportedMolecule) else None
        for molecule in result.molecules
    )
    imported_metadata = tuple(value for value in metadata if value is not None)
    shared_metadata = (
        imported_metadata[0]
        if imported_metadata
        and len(imported_metadata) == len(metadata)
        and all(
            value.structural_input == imported_metadata[0].structural_input
            and value.conformers == imported_metadata[0].conformers
            for value in imported_metadata
        )
        else None
    )
    requested = dict(result._requested_payload)
    requested["structural_input"] = (
        None if shared_metadata is None else shared_metadata.structural_input
    )
    requested["conformers"] = None if shared_metadata is None else shared_metadata.conformers
    return _native_adapters._dumps(
        result._native,
        result.molecules._native_molecules,
        identities,
        tuple(() if value is None else value.diagnostics for value in metadata),
        requested,
        format,
        sdf_version or "v3000",
    )


def write(
    path: str | PathLike[str],
    result: CalculationResult,
    *,
    format: OutputFormat,
    sdf_version: SdfVersion | None = None,
) -> None:
    """Serialize a calculation result to a UTF-8 text file."""

    if not isinstance(path, (str, PathLike)):
        raise TypeError("path must be a string or path-like value")
    Path(path).write_text(
        dumps(
            result,
            format=format,
            sdf_version=sdf_version,
        ),
        encoding="utf-8",
    )


__all__ = [
    "InputFormat",
    "OutputFormat",
    "INPUT_FORMATS",
    "OUTPUT_FORMATS",
    "SdfVersion",
    "RecordSelection",
    "BondStrategy",
    "ConformerSelection",
    "parse",
    "read",
    "dumps",
    "write",
]
