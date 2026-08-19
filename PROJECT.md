# ChargeFW Project Guide

This document is the source of truth for ChargeFW's technical state, compatibility context, and
product roadmap. Read [AGENTS.md](AGENTS.md) for implementation rules, [TODO.md](TODO.md) for
actionable deliverables, and [README.md](README.md) for build and test commands.

## Purpose

ChargeFW is a C++23 framework for empirical partial atomic-charge calculation. It provides a
toolkit-neutral core model, built-in empirical methods, parameter-set loading and classification,
method applicability checks, and structured charge results.

It is being developed as a modern, library-first successor to **ChargeFW2**, the current
computational engine of [Atomic Charge Calculator III (ACC III)](https://acc.biodata.ceitec.cz).
ACC III is described in Raček *et al.*, *Atomic Charge Calculator III: a modern platform for
calculating partial atomic charges*, **Nucleic Acids Research** (2026),
DOI [10.1093/nar/gkag379](https://doi.org/10.1093/nar/gkag379).

ChargeFW is not yet the ACC III backend. `old/` is an archived ChargeFW2 copy used for behavior,
parameter, and compatibility research.

## Current architecture

```text
External inputs and front ends                    Only the demo CLI exists today
    (current: native molecular I/O and Gemmi; planned: Python/NumPy, RDKit, JSON/WASM)
                              |
                              v
core::MoleculeCollection
    |- core::Molecule: atoms, bonds, conformers
    |- core::Atom: atomic number, formal charge, source name
    |- core::Bond: atom indices and order
    `- core::Conformer: coordinates sharing the molecule topology
                              |
                              v
features::PreparedMoleculeCollection
    |- PreparedMolecule: molecule + cached TopologyFeatures
    `- ConformerFeatures: on-demand geometry view per conformer
                              |
                 +------------+------------+
                 |                         |
                 v                         v
     parameters::ParameterSet       methods::MethodRegistry
                 |                         |
                 v                         v
      ParameterClassification     MethodRequirements/options
                 |                         |
                 +------------+------------+
                              v
                 methods::find_applicable_methods()
                              |
                              v
                      methods::ApplicableMethod
                              |
                              v
             calculation::calculate(CalculationRequest)
                               |
                               v
                       charges::ChargeSet

Normal application path: calculation::calculate() composes applicability,
deterministic candidate selection, and calculate_charges().
```

### Important public types

| Area | Types | Role |
|---|---|---|
| Core | `Atom`, `Bond`, `Conformer`, `Molecule`, `MoleculeCollection` | Input molecular graph and coordinates. |
| Features | `TopologyFeatures`, `ConformerFeatures`, `PreparedMolecule` | Cached/derived topology and geometry. |
| Parameters | `ParameterSet`, `ParameterClassification`, `ParameterView` | Parameter storage, matching, and method-facing lookup. |
| Methods | `Method`, `MethodRegistry`, `MethodRequirements`, `MethodOptions`, `ApplicableMethod` | Algorithm interface, capabilities, selection, and execution. |
| Charges | `AtomicCharges`, `ChargeAssignment`, `ChargeSet`, `ChargeCollection` | Atom-indexed calculated results and provenance. |
| Calculation | `CalculationRequest`, `CalculationResult`, `ApplicationCalculationRequest` | Selected-candidate execution and application-facing automatic facade. |

## Explicit native workflow

```cpp
core::MoleculeCollection molecules{/* validated molecules */};
features::PreparedMoleculeCollection prepared{molecules};

const auto parameter_sets = parameters::load_default_parameter_sets();
const auto& registry = methods::method_registry();

std::vector<const methods::Method*> candidates;
for (const auto& method : registry.methods()) {
    candidates.push_back(method.get());
}

const auto applicability =
    methods::find_applicable_methods(
        {.molecules = prepared,
         .methods = candidates,
         .parameter_sets = parameter_sets,
         .classification_options = {.permissive_types = false}});

// Application policy may select a reported candidate itself, or use the deterministic helper.
const auto* selected = calculation::select_applicable_method(applicability);
if (selected == nullptr) {
    // Inspect applicability.rejected.
} else {
    const auto result = calculation::calculate({.molecules = prepared, .selected = *selected});
}
```

The explicit path is:

```text
PreparedMoleculeCollection
  + ApplicabilityRequest { methods, parameter sets, classification options }
  -> ApplicabilityResult { applicable candidates with resolved classifications, rejected candidates }
  -> caller selection or select_applicable_method()
  -> CalculationRequest { prepared molecules, selected ApplicableMethod }
  -> CalculationResult { ChargeSet }
```

Classification options belong to applicability. Once a candidate exists, calculation constructs
`ParameterView` from that candidate's stored atom/bond classifications and never reclassifies.

## Convenience application workflow

For normal application use, `calculation::calculate(ApplicationCalculationRequest)` provides the
owned convenience facade:

```text
ApplicationCalculationRequest {
  molecules,
  parameter sets,
  optional method ID,
  optional parameter-set ID,
  classification options
}
  -> prepare
  -> applicability
  -> deterministic selection
  -> execution
  -> ApplicationCalculationResult { optional ChargeSet, applicability diagnostics }
```

The facade consumes `classification_options` only while determining applicability. It selects the
applicable candidate with the highest method priority, then the highest parameter-set priority. Ties
are resolved deterministically by method ID and parameter-set ID. Higher priorities therefore denote
maintainer-curated automatic preference, not a universal scientific quality ranking. Explicit IDs
restrict the candidates and fail if unavailable or inapplicable rather than silently falling back.
The result retains applicability diagnostics when no candidate can be calculated.

The current `chargefw` executable is a molecular-file demonstration. It autodetects `.sdf`, `.mol`,
`.mol2`, `.pdb`, `.cif`, `.mmcif`, and ChargeFW `.json` input from the file extension, rejects the
entire input on the first malformed record, loads bundled parameter sets, and autodetects the
highest-priority applicable method and parameter set. All inputs produce JSON and mmCIF outputs.
Native-molecular and JSON input also produce SDF and MOL2; structural PDB/mmCIF input does not.
Same-format SDF and MOL2 outputs
preserve the source; other molecular outputs are generated. The required output-directory argument is
created when absent, and output filenames use `<input-stem>.chargefw`. JSON molecules with multiple
conformers are rejected because one record cannot currently represent all assignments consistently.
For PDB/mmCIF input, `--structural-selection` selects all records, polymers plus ligands excluding
water, or polymers only; `--structural-bonds` selects no bonds, explicit connectivity, compact
templates, or their hybrid. Both options reject non-structural formats. It is not yet a full
user-facing file/SMILES CLI.

### Calculation granularity

`Method::calculate()` operates on one molecule and zero or one conformer. The collection-level
`methods::calculate_charges()` operation applies one selected method/parameter candidate across the
prepared collection: it calculates every conformer of every molecule for geometry-dependent
methods, and one conformer-independent assignment per molecule otherwise. A molecule without a
conformer remains valid for methods that do not require geometry; geometry-dependent methods report
it as inapplicable. Future streaming and batch helpers must preserve molecule, conformer, and atom
order and report each target independently.

### Calculation selection and result contract

`methods::find_applicable_methods(ApplicabilityRequest)` evaluates supplied method and parameter-set
candidates and returns resolved `ApplicableMethod` candidates. Each candidate stores its atom/bond
parameter classifications. `calculation::select_applicable_method()` applies the deterministic
priority ordering, or callers can select a candidate themselves. `calculation::calculate(
CalculationRequest)` executes that selected candidate using its stored classifications; it does not
accept or reevaluate classification policy. A successful `CalculationResult` contains one
`charges::ChargeSet`; its owned method ID and optional parameter-set ID identify the candidate used.

`ApplicationCalculationResult::applicability` retains considered applicable and rejected candidates
for the convenience facade. If no candidate is applicable, its `charges` is empty; calculation
failures after selection are reported as failures rather than silently treated as inapplicability.

`ChargeSet` preserves calculation targets as `ChargeAssignment` entries. For geometry-dependent
methods, it contains one assignment for every conformer of every input molecule, each identified by
its molecule index and conformer index. For geometry-independent methods, it contains one
assignment per molecule and no conformer index. Each assignment's atomic-charge vector follows the
source molecule's atom order.

`CalculationRequest` is intentionally a low-level native execution view: it contains a prepared
collection and a selected `ApplicableMethod`, both owned by the caller. The selected candidate must
come from applicability; its stored classifications are the only parameter mapping used for
execution. `ApplicationCalculationRequest` is the binding-friendly owned facade: it accepts a
native molecule collection and parameter sets, resolves registered methods, applies classification
policy, and selects supplied candidates by optional IDs. Omitted IDs use deterministic automatic
selection; explicit unavailable or inapplicable IDs fail rather than silently falling back.
Method-specific options remain a low-level advanced-native feature until their application-facing
policy is specified. All integrations must compose these calculation paths rather than reimplement
selection or scientific behavior.

### Scalable execution policy

Large-molecule acceleration is a primary product capability, not deferred work. The application
calculation contract should expose an execution policy distinct from scientific method options:
`full` performs the exact/reference calculation, while `cutoff(radius)` and `cover(radius)` perform
explicitly reported fragment approximations. Keeping execution policy outside individual `Method`
implementations allows the same spatial decomposition, fragment mapping, parallel scheduling, and
charge-reconciliation machinery to be reused by EEM/QEq-like and SQE-family methods.

The public `ExecutionPolicy`, `ExecutionSelection`, and `ResourcePolicy` value types establish the
validated vocabulary: concrete `full`, `cutoff(radius)`, and `cover(radius)` policies; an `automatic`
caller preference; and a shared full-atom threshold with an unlimited representation. Reduced radii
must be finite and at least 8 angstrom; full rejects a radius. `ApplicationCalculationRequest` now
carries an execution selection, but it is preparatory only: until execution availability and concrete
plan selection are implemented, the facade retains its existing full-calculation behavior and does
not claim cutoff or cover support.

ChargeFW2 implemented cutoff and cover only through `EEMethod`. It also switched modes silently at
fixed atom-count thresholds. ChargeFW should first reproduce and validate that behavior for the
methods that previously supported it, but must not retain silent switching. Any future `auto`
execution policy needs benchmark-derived thresholds and must report the selected mode and radius.

Extending cutoff and cover to SQE, SQE+q0, and SQE+qp is an intentional new approximation rather
than a compatibility claim. Fragment construction must preserve induced topology, bond parameters,
atom mapping, and each method's initial-charge semantics. Fragment target-charge allocation,
overlap reconciliation, and final total-charge correction must be specified centrally and tested
against full calculations as radius increases. Methods must declare whether they support each
execution policy; unsupported combinations fail explicitly.

Spatial neighbor search and reusable fragment data belong in `features`, not `core::Molecule`.
Validation must cover charge conservation, deterministic results, atom-order preservation,
convergence toward `full`, accuracy/error envelopes, runtime, memory, and parallel execution.

## Integration and distribution architecture

Python is expected to be the primary workflow integration, while the C++ library and C++ CLI remain
supported products. All front ends must converge on the same native request/result behavior:

```text
RDKit Mol ── pure-Python converter ─┐
Biopython/NumPy ─ Python converter ─┼─> toolkit-neutral molecule data
SDF/Gemmi/C++ caller ─ C++ adapter ─┘             |
                                                   v
                                      native calculation facade
                                                   |
                                                   v
                                  charges + mapping + provenance
```

### Python and RDKit

The first RDKit integration should be implemented at the Python level. A Python converter can read
an ordinary `rdkit.Chem.Mol` and pass atomic numbers, formal charges, indexed bonds, and conformer
coordinates through the native binding. This gives users a direct `calculate(mol)` workflow while
avoiding RDKit C++ headers, link libraries, ABI compatibility, and wheel-loading concerns. The
conversion cost is expected to be small compared with classification and numerical calculation and
must be measured before introducing a compiled RDKit dependency.

The base Python wheel should contain the ChargeFW native extension, toolkit-neutral Python/NumPy
API, and bundled parameter data. `uv add chargefw` must work without RDKit;
`uv add "chargefw[rdkit]"` should add the optional Python RDKit dependency and convenience module.
RDKit must be imported lazily so base installations remain functional. A direct C++ adapter for
`RDKit::ROMol` remains a valid future optional target, but only for a demonstrated native consumer
and never as a dependency of `chargefw_core` or the base wheel.

RDKit import is a representation conversion, not an implicit chemistry-preparation step. It must
preserve source atom order, formal charges, bond information, and conformer identity. Sanitization,
hydrogen addition/removal, protonation, tautomer changes, embedding, and optimization require
separate opt-in helpers and explicit provenance. Writing partial charges back to an RDKit molecule
is also a separate, explicit mutation operation; partial charges must never replace formal charges.

### Other integrations

NumPy-style arrays form the toolkit-neutral Python interchange and make custom scientific, ML, and
simulation workflows possible without package-specific bindings. Biopython remains a possible
secondary adapter for structural-biology object models, but coordinates and hierarchy do not
guarantee a complete chemical graph or bond orders; it must map onto explicit connectivity, model,
component, and alternate-location policies. Gemmi 0.7.4 currently backs native PDB and mmCIF readers in
`chargefw_core`. The PDB reader maps compatible models to conformers of one molecule; the mmCIF
reader parses the CIF document eagerly, then converts each coordinate-bearing `data_` block lazily
to a separate record, with compatible models in that block represented as conformers. Both readers
preserve selected atom order, atom names, formal charges, conformer identity, and record identity.
They select blank alternate locations before `A` and then the first occurrence, and support
all-records, polymers-and-ligands (water excluded), and polymers (HETATM excluded) selection modes.

Connectivity is an explicit input option with four strategies: `none` (the default), `templates`,
`explicit_bonds`, and `hybrid`. The compact built-in CCD-derived template catalog covers standard
amino acids, standard RNA/DNA nucleotide names, and water, and adds sequential peptide C-N and
nucleotide O3'-P links. PDB explicit connectivity reads `CONECT` plus covalent/disulfide Gemmi
connections; mmCIF explicit connectivity reads local `_chem_comp_bond` rows plus covalent/disulfide
`_struct_conn` connections. `hybrid` combines explicit and template connectivity, deduplicating atom
pairs and retaining an explicit bond order on conflicts. Distance-based perception and an external
full-CCD provider are not implemented.

An exploratory full-CCD packed-template benchmark represented 50,782 components and 2,523,648
bonds as pooled atom names, a flat bond table, and a sorted component index. The optimized stripped
lookup executable was about 1 MB, but compiling the generated 41 MB C++ source was demanding. This
keeps a separately loadable full-CCD provider viable while compact built-in templates remain the
normal path.

Open Babel, MDAnalysis, MDTraj, and OpenFF may receive thin convenience bridges when concrete
workflows justify them. They should normally convert through the toolkit-neutral Python boundary
rather than becoming native core dependencies. The confirmed standalone CLI and bounded-memory
batch workflow justify a focused native SDF adapter; broad format support should otherwise come from
optional established toolkits.

### Molecular I/O

An application requires a real-molecule input path even when no external toolkit is installed. The
standalone C++ CLI should therefore gain a bounded-memory native MOL/SDF stream adapter as its first
file interface, beginning with a clearly documented subset and explicit rejection of unsupported
V2000/V3000 features. It must preserve record identity, atom order, formal charges, bonds, and
coordinates and return record-scoped errors without silently repairing chemistry. SDF output should
retain the source mapping and attach charges and provenance without changing molecular semantics.

The native SDF path is not intended to become a universal chemistry toolkit. Python RDKit provides
broad SMILES/SDF/Mol2 workflows; the required Gemmi adapter and future optional Biopython adapters
cover structural biology. Each adapter translates to the same toolkit-neutral molecule and
calculation contracts, so adding a format or package cannot introduce a separate selection or
scientific-policy implementation.

#### Adapter record contract

Adapters exchange `adapters::ImportedMoleculeRecord` rather than attaching format state to
`core::Molecule`. Each successful record owns its native molecule, source identity, source-to-native
atom/conformer mapping, and non-fatal diagnostics. A mapping is explicit even when it is the identity
mapping, so Gemmi and future RDKit and Python adapters can preserve source order without adapter-
specific result rules. An adapter may retain an opaque, format-tagged source payload when its writer
must enrich an existing record rather than reconstruct it. The future Gemmi mmCIF writer will use
this to preserve structural categories and append the `_sb_ncbr_partial_atomic_charges_meta` and
`_sb_ncbr_partial_atomic_charges` loops from the archived exporter. Record failures use
`adapters::MoleculeRecordError`, which retains identity, message, and optional source line for
reporting the first malformed record before terminating import. This deliberately defines no common
file handle, property bag, hierarchy model, or chemistry repair policy: those remain adapter-specific
until a concrete format requires them.

Native adapters are named by format and direction. The current `json_input`, `mol_input`,
`sdf_input`, `mol2_input`, and Gemmi-backed `pdb_input` adapters import molecule records;
`mmcif_input` parses the Gemmi CIF document eagerly and lazily converts each coordinate-bearing
mmCIF `data_` block as a record; `json_output::JsonWriter` serializes the format-neutral
`ChargeResultDocument` used by the CLI and future integrations.
`mol2_output::Mol2Writer` preserves a source MOL2 file while replacing or adding atom partial-charge
fields for one selected record, or generates a basic MOL2 record from native graph and conformer
data when the input is another format. Generated MOL2 uses element-symbol atom types and a single
`CHARGEFW` substructure; it does not infer Tripos typing or source-specific substructure semantics.
`sdf_output::SdfWriter` preserves source SDF records and writes atom-order charge vectors as numbered
`CHARGEFW_CHARGES_<type-id>` data fields. Replace mode removes existing fields with that owned prefix
before writing the new set; append mode retains existing fields and adds another set. It can also
generate an explicit V2000 or V3000 SDF record from native graph and conformer data, preserving
formal charges in the MOL representation and partial charges in SDF properties.
Additional writers must consume that export model rather than implementing calculation-result
serialization in a front end.

#### Format and connectivity policy

Test molecular files are organized under `tests/fixtures/`: small hand-authored cases belong in
`synthetic/<format>/`, while intact real-world inputs belong in `corpus/<format>/<subject>/` with a
short provenance note. Future PDB and mmCIF fixtures should use this same format-first structure,
with subdirectories for scenarios such as multi-model, multi-component, and alternate locations.

The native adapter scope is intentionally narrow: MOL/SDF, Tripos MOL2, and a versioned ChargeFW JSON
document are dependency-free CLI input paths; Gemmi-backed PDB/mmCIF readers are library APIs but are
not yet wired into the demo CLI. JSON documents use `schema_version: "1.0"` and a
`molecules` array. Each molecule has optional `id` and `name`, required `atoms` entries with
`atomic_number` and `formal_charge`, optional indexed `bonds`, and optional conformers with coordinate
triplets; atom names are intentionally excluded. Array order is authoritative for atom and conformer
mapping, and connectivity is never inferred from coordinates. The initial MOL/SDF reader supports V2000 atom and bond blocks plus `M  CHG`
formal-charge records, and
V3000 CTAB `COUNTS`, `ATOM`, and `BOND` blocks with `CHG=` attributes. It imports aromatic bond order
as single bonds
and records ignored `CFG=` stereochemical attributes as diagnostics; it rejects query atoms, unknown
elements, unsupported bond orders, unsupported V2000 properties, and unsupported V3000 attributes.
`parse_mol()` handles one standalone MOL record through `M  END`; `SdfReader` adds bounded-memory
multi-record framing and skips SDF data fields without interpreting them. `Mol2Reader` supports the
MOL2 `MOLECULE`, `ATOM`, and `BOND` sections, standard element-prefixed atom types, and numeric bond
types; aromatic bond types are imported as single bonds. MOL2 partial-charge fields are ignored, never treated as formal charges, and
reported once per record when nonzero values are present. It is not a general chemistry toolkit. RDKit
support is optional and normally
enters through the Python bridge; if a native RDKit adapter is later justified, it is a separately
selected backend and never silently replaces the native SDF reader. PDB/mmCIF parsing uses
Gemmi-backed adapters; mmCIF charge export is also Gemmi-backed. The project will not
implement those formats itself.
`pdb_input` treats compatible PDB models as conformers of one molecule, retaining atom names,
formal charges, and coordinates. It selects blank alternate locations before `A`, then the first
occurrence, supports all-records, polymers-and-ligands (water excluded), and polymers
(HETATM excluded) selection modes. `mmcif_input` applies the same selection and alternate-location
rules, represents compatible models inside one `data_` block as conformers, and returns separate
records for separate coordinate-bearing `data_` blocks. Both readers expose `none`, `templates`,
`explicit_bonds`, and `hybrid` connectivity strategies as described above.
The Gemmi writer semantically preserves the parsed mmCIF document and appends or replaces the
SB-NCBR charge dictionary/categories without reconstructing unrelated structural data. PDB input is
converted through Gemmi. Other inputs generate one self-contained `UNL` component block per molecule
record. Gemmi serialization may normalize presentation and is not byte-for-byte preservation.

#### Preservation-oriented output

When output format matches a filesystem-backed input, writers should transform a byte-for-byte copy
of the source file rather than reconstructing a molecule document. Preservation means making only
the smallest required format-specific edits: SDF writers append ChargeFW-owned data fields before
the record delimiter; mmCIF writers append ChargeFW categories inside the selected `data_` block;
MOL2 writers replace the existing ATOM partial-charge token or add the missing optional trailing
fields required to attach a charge; and future PQR writers replace mapped atom charge values while
retaining radii. A writer rejects only when it cannot safely locate or validate the selected source
record/block and mapped atoms, never merely because an optional charge field is absent. Generated
output is the explicit alternative for an output format that differs from input or lacks a readable
source file.

Structural adapters must select connectivity explicitly. The implemented policies are `none`,
`templates`, `explicit_bonds`, and `hybrid`; the latter combines compact built-in component templates
for intra-component and sequential polymer bonds with explicit PDB/mmCIF records. A future
distance-perception policy may infer connectivity from covalent radii and coordinates, but it must be
opt-in, never a fallback, and must report inferred bond orders and provenance. Alternate locations
likewise use an explicit deterministic policy: blank alternate location, then `A`, then the first
occurrence. Adapters preserve selected native ordering; omitted-altloc mapping and diagnostics remain
to be completed.

A method may not require topology directly, but parameter classification can require atomic
environment or highest bond order. Import success therefore does not imply calculation applicability:
the existing method/parameter prerequisite checks remain responsible for rejecting insufficient or
unsupported topology. Every structural calculation/export result must retain the selected
connectivity and alternate-location policies and their provenance.

### Packaging boundaries

Native CMake installation and Python distribution are independent delivery paths over the same C++
implementation. Native consumers should eventually use exported CMake targets and installed
parameter data. Python users should receive nanobind-based binary wheels built from `pyproject.toml`
with a CMake-aware backend such as scikit-build-core, containing the extension and parameter
resources; they must not need a prior system ChargeFW installation or `CHARGEFW_PARAMETER_DIR`.
Conda-forge is a complementary channel, especially for environments already using RDKit, not a
replacement for ordinary `uv`/pip wheels. A versioned OCI/Docker image should provide the standalone
CLI and bundled data for reproducible batch deployments.

JSON is the intended stable serializable and process boundary for the CLI, WebAssembly, services,
and MCP, but it is not required merely to pass an in-process RDKit molecule to the native extension.
Both the Python object API and JSON serialization must expose the same calculation targets, source
atom mapping, selected method and parameters, options, warnings, and diagnostics.

## Implemented methods and parameters

The current registry contains 22 methods:

```text
abeem, charge2, delre, denr, dummy, eem, eqeq, eqeqc, formal, gdac,
kcm, mgc, mpeoe, peoe, qeq, sfkeem, smpqeq, sqe, sqeq0, sqeqp, tsef, veem
```

Bundled JSON parameter sets cover these parameterized methods and variants. `data/parameters/`
is installed under `share/chargefw/parameters`.

### Compatibility gap with ChargeFW2

The archived ChargeFW2 registry contains the same 22 methods. All nine archived SQE-family
parameter sets have been migrated. Numerical parity remains a release-blocking objective because
ACC III highlights SQE+qp.

## ChargeFW2 research summary

ChargeFW2 is a functional, application-oriented engine with CLI/Python bindings, custom
SDF/Mol2 reading, Gemmi PDB/mmCIF reading, output writers, OpenMP, and a nanoflann KD-tree. It
has useful production behavior, but combines concerns that are intentionally separated here:

| ChargeFW2 pattern | ChargeFW replacement |
|---|---|
| Atom stores graph data, coordinates, PDB residue data, Mol2 types, and classification state. | Core graph stays minimal; adapters and derived feature layers own source-specific metadata/caches. |
| Molecule owns all-pairs topology matrices and spatial KD-tree. | `TopologyFeatures` and `ConformerFeatures` own derived data. |
| Parameter classification mutates atoms and bonds. | Immutable `ParameterClassification` and `ParameterView`. |
| Method owns mutable parameters pointer and option values. | Stateless `Method::calculate(CalculationInput)`. |
| Suitability is mainly filtering/console behavior. | Structured applicability and prerequisite diagnostics. |
| PDB/mmCIF reader uses the first model only. | Core can represent multiple conformers; adapters should define model/altloc policy. |
| `full`, `cutoff`, and `cover` are coupled inside `EEMethod` with silent size-based switching. | Reusable execution policy is explicit, capability-checked, benchmarked, and provenance-rich. |

The new implementation must preserve validated scientific behavior before making intentional
improvements. Compatibility fixtures should use identical topology, coordinates, formal charge,
options, and parameter data.

### SQE-family research

ChargeFW2 implements the SQE family with a signed bond-incidence matrix `T`. SQE solves
`(T A Tᵀ + diag(kappa)) p = T b`, then returns `q = Tᵀ p`, where `A` contains atom hardness on the
diagonal and the Gaussian-width Coulomb interaction `erf(d / sqrt(2 wi² + 2 wj²)) / d` off the
diagonal, and `b = -electronegativity`. SQE therefore conserves zero total charge within each
connected component. `sqeq0` adds formal charges as initial charges; `sqeqp` uses parameterized
`q0` after uniformly correcting its total to the molecular formal charge. All require coordinates,
atom parameters `electronegativity`, `hardness`, and `width`, plus bond parameter `kappa`; SQE
does not itself use formal charges, while SQE+q0 requires them. SQE+qp additionally requires atom
parameter `q0`, which it uniformly corrects to the molecule's formal-charge total. The archived
implementation does not explicitly diagnose singular or ill-conditioned transfer systems.

## Product direction

The product work is organized around four outcomes: scalable `full`/`cutoff`/`cover` calculation;
real-molecule I/O and a useful C++ CLI; modern adapters led by Python/NumPy and RDKit; and easy,
reproducible installation through wheels, Conda, native CMake packages, and a container image. The
native request/result/provenance contract is the shared foundation for all four. Scientific
comparison, approximation validation, diagnostics, and release automation are cross-cutting quality
gates rather than separate product directions.

The prioritized delivery checklist and validation criteria are maintained in [TODO.md](TODO.md).
WebAssembly, MCP, ML/GNN adapters, and additional toolkit bridges must not be pre-built ahead of
validated need.

## Scientific and integration principles

- Preserve source atom indexing and map every output charge back to its source atom.
- Record method, parameter set, method options, conformer, approximation policy, and warnings in
  every serializable result.
- Keep exact/reference calculations distinct from cutoff/cover approximations.
- Make automatic selection explainable and overridable.
- Treat empirical charge results as method-dependent estimates; expose coverage and numerical
  limitations rather than implying universal validity.
