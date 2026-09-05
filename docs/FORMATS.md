# Molecular input and charge output formats

This document describes the format behavior shared by ChargeFW's native adapters, command-line
application, and language bindings. API names and command-line file-selection policy are documented in
[Native C++ library](NATIVE.md), [Python package](PYTHON.md), and
[Command-line interface](CLI.md), respectively.

ChargeFW imports files into a toolkit-neutral molecule containing source-ordered atoms, bonds, formal
charges, and zero or more Cartesian conformers. Import does not sanitize molecules, infer arbitrary
bonds, change protonation or hydrogen count, generate conformers, or optimize geometry. Atom and molecule
order is retained in charge assignments and in every output mapping.

## Format overview

| Format | Input | Charge output | Record model |
| --- | --- | --- | --- |
| MOL V2000/V3000 | Yes | Through generated SDF | One molecule, one conformer |
| SDF V2000/V3000 | Yes | Preserved or generated SDF | One molecule per record, one conformer each |
| Tripos MOL2 | Yes | Preserved or generated MOL2 | One molecule per `MOLECULE` record, one conformer each |
| ChargeFW molecule JSON 1.0 | Yes | No; result JSON uses a different schema | `molecules` array, zero or more conformers each |
| PDB | Yes, through Gemmi | Converted to mmCIF | One molecule, models become conformers |
| mmCIF | Yes, through Gemmi | Preserved or generated mmCIF | One molecule per coordinate-bearing block |
| ChargeFW result JSON 1.0 | No | Yes | One result record per imported molecule |

Readers report malformed or unsupported molecular data as errors. A successful imported record also
carries its source name, zero-based record index, format-derived record ID, and non-fatal diagnostics.
Streaming SDF and MOL2 readers consume one record at a time; molecule JSON, PDB, and mmCIF parsing retains
the complete source document in memory.

## MOL and SDF input

The MOL reader accepts one V2000 or V3000 record. The SDF reader applies the same MOL parser to each
record and skips SDF data fields after `M  END`; data fields do not become molecule properties.

The supported subsets import:

- concrete element symbols represented by ChargeFW's periodic table;
- Cartesian coordinates as one conformer;
- source atom order and single, double, or triple covalent bonds; and
- formal charges from V2000 `M  CHG` records or V3000 atom `CHG=` attributes.

Query and wildcard atoms are rejected. Bond order 4 is imported as a single bond. V3000 coordinate bonds
(order 9) and hydrogen bonds (order 10) are omitted because the core graph cannot represent them; one
warning per omitted order is attached to the molecule record. Other trailing V3000 atom and bond
attributes are accepted but not used.

V2000 properties other than `M  CHG` are not interpreted. One warning per ignored three-character
property code is attached to the molecule record. Repeated ignored-property and omitted-bond warnings are
coalesced within a record and retain the first record-relative line number.

The first MOL header line supplies the record ID and, unless another name is available, the molecule
name. MOL and SDF input always produces one conformer named `input`.

## MOL2 input

The MOL2 reader consumes `@<TRIPOS>MOLECULE` records and reads their `MOLECULE`, `ATOM`, and `BOND` data.
Unrelated sections are skipped before the atom section. The declared atom and bond counts, positive unique
atom IDs, bond references, coordinates, and required sections are validated.

The element is taken from the prefix before `.` in a standard atom type, while the atom-name column is
retained. Dummy, lone-pair, wildcard, halogen-group, hetero-group, and heavy-atom-group types are rejected.
Numeric bond types 1, 2, and 3 are imported directly; numeric type 4 and aromatic type `ar` become single
bonds.

MOL2 atom charges are partial charges, not formal charges. They are parsed for validity but are not used;
all imported formal charges are zero. A record containing any nonzero input partial charge receives a
`partial_charges_ignored` warning. Each record produces one conformer named `input`.

## ChargeFW molecule JSON input 1.0

Molecule input JSON is an explicit interchange format and is distinct from ChargeFW result JSON. A
document has this shape:

```json
{
  "schema_version": "1.0",
  "molecules": [
    {
      "id": "water-1",
      "name": "water",
      "atoms": [
        {"atomic_number": 8, "formal_charge": 0},
        {"atomic_number": 1, "formal_charge": 0}
      ],
      "bonds": [
        {"atoms": [0, 1], "order": 1}
      ],
      "conformers": [
        {"id": "model-1", "coordinates": [[0, 0, 0], [0.96, 0, 0]]}
      ]
    }
  ]
}
```

`schema_version`, `molecules`, and each molecule's non-empty `atoms` array are required. Every atom
requires integer `atomic_number` and `formal_charge` members. Optional `bonds` use zero-based atom indices
and bond orders 1, 2, or 3. Optional `conformers` contain one finite three-number coordinate row per atom;
their string `id` becomes the conformer name. Molecule `id` and `name` are optional strings. Unknown
members are ignored, but atom names are intentionally not part of schema 1.0.

Array order defines molecule, atom, bond, and conformer order. Readers may select either the first
conformer or all conformers. The whole JSON document is parsed before records are returned; malformed
JSON, another schema version, incorrect types, invalid indices, and inconsistent molecule dimensions are
errors.

## PDB and mmCIF input

PDB and mmCIF are parsed through Gemmi. PDB produces one molecule record. Each mmCIF data block containing
`_atom_site.id` produces one molecule record; blocks without coordinate data are skipped. The source block
name is retained as the mmCIF record ID.

The first selected model defines atom order, elements, formal charges, names, topology, and source
mapping. Later models become conformers only when they have the same selected atom count, element, formal
charge, and atom-name sequence. Readers may retain the first model or all models. Unknown elements, empty
selections, absent models, and incompatible model sequences are errors.

### Structural selection

Selection is applied before conformer validation and bonding:

- `all` retains polymer, ligand, and water residues;
- `polymers-and-ligands` excludes `HOH` water residues; and
- `polymers` excludes all hetero residues.

Only one alternate location is imported for each atom name in a residue. A blank location is preferred,
then `A`, then the first location present in source order.

### Connectivity

Structural readers provide four explicit bond strategies:

- `none` imports no bonds;
- `explicit` imports PDB `CONECT` and covalent/disulfide connections, or mmCIF component bonds and
  covalent/disulfide connections;
- `templates` applies ChargeFW's built-in residue templates and peptide/nucleotide polymer links; and
- `hybrid` combines explicit and template bonds, with explicit connectivity taking precedence.

Explicit structural aromatic bonds become single bonds. Templates cover common amino acids,
nucleotides, water, and basic sequential peptide and nucleotide links; they are not a complete Chemical
Component Dictionary. No distance-based bond perception is performed.

## Charge output

All serialized partial charges are rounded to four decimal places. Internal native and Python result
objects retain calculation precision. Writers validate assignment order, cardinality, and source
references where the output representation provides them rather than silently reordering assignments.

### SDF

Preservation-oriented output copies every source record and its unrelated fields. In the default
`replace` mode, existing numbered `CHARGEFW_CHARGES_*` and `CHARGEFW_CHARGE_METADATA_*` fields are removed;
`append` mode retains them. Each charge set is written as:

```text
> <CHARGEFW_CHARGES_1>
-0.8000 0.4000 0.4000

> <CHARGEFW_CHARGE_METADATA_1>
type=empirical; method=eem; parameter_set=example; software_name=ChargeFW; software_version=...
```

Generated output can use V2000 or V3000. It writes the native atom order, formal charges, supported bond
orders, and the conformer identified by the charge assignment. V2000 generation is limited to 999 atoms
and 999 bonds. Generation requires coordinates and one assignment per numbered charge property; multiple
properties must target the same conformer.

### MOL2

Preservation-oriented output retains source text and replaces the ninth atom field with each calculated
charge. If source atom rows have no substructure and charge fields, `1 UNL` and the charge are added.

Generated output writes a `SMALL` molecule with `USER_CHARGES`, native atom and bond order, one selected
conformer, and one charge assignment. Generated atom types are element symbols only; ChargeFW does not
infer Tripos atom types or substructures. Both preserved and generated output require exactly one
assignment for each molecule.

### mmCIF

mmCIF charge output follows the SB NCBR partial atomic charges dictionary version 1.1. It adds
`_sb_ncbr_partial_atomic_charges_meta` metadata and `_sb_ncbr_partial_atomic_charges` values, plus the
dictionary declaration in `_audit_conform`. Charge rows refer to `_atom_site.id`, preserving the mapping
of structural selections and conformers.

For mmCIF input, preservation-oriented output retains the source document and unrelated categories. Its
default `replace` mode replaces existing ChargeFW charge categories; `append` allocates subsequent numeric
charge-type IDs. PDB input is converted to mmCIF after applying the imported structural selection.
Unselected PDB alternate locations are omitted by that conversion, while an mmCIF source document retains
unselected rows and attaches charges only to selected atom IDs.

Generated mmCIF creates one non-polymer `UNL` data block per molecule, including elements, formal charges,
single/double/triple bonds, and all native conformers. Block and atom IDs are made unique when necessary.
Conformer-independent assignments are attached to every retained conformer; conformer-specific
assignments are attached only to their corresponding model. Coordinates are required.

### ChargeFW result JSON 1.0

Result JSON is the complete machine-readable calculation record. Its top level contains:

- `schema_version`, generator identity, overall status, and document diagnostics;
- source-ordered `results`, each with source identity, record status, diagnostics, and successful charge
  assignments; and
- optional `calculation_provenance` containing requested and effective calculation policy and execution
  metrics.

Each assignment declares `atom_order: "source"`, contains a source-ordered `charges` array and its
serialized `total_charge`, and includes a zero-based `conformer_index` when the calculation is
conformer-specific. Failed or cancelled records omit assignments. Diagnostics have stable severity, code,
and message fields and may include zero-based molecule, atom, bond, or conformer indices and a one-based
source line number.

Requested provenance records method and parameter selection, classification mode, method options,
available conformer and structural import policy, execution request, resource thresholds, and thread
limit.
Effective provenance records the resolved method, parameter set, complete options, execution policy, and
warnings. When supplied by the application, metrics include UTC start/end timestamps, parsing,
applicability, computation, writing and total runtimes, and peak resident memory. Durations and memory are
rounded to three decimal places.
