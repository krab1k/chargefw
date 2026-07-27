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
    (planned: Python/NumPy, Python RDKit, SDF, Gemmi, JSON/WASM)
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
                   methods::calculate_charges()
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
| Calculation | `CalculationRequest`, `CalculationResult` | High-level applicability, deterministic selection, and execution facade. |

## Typical library workflow

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
    methods::find_applicable_methods(prepared, candidates, parameter_sets);

// Application policy selects an explicitly reported candidate.
const auto& selected = applicability.applicable.front();
const charges::ChargeSet result = methods::calculate_charges(selected, prepared);
```

For normal application use, `calculation::calculate()` provides autodetection over the supplied
method and parameter candidates. It selects the applicable candidate with the highest method
priority, then the highest parameter-set priority. Ties are resolved deterministically by method ID
and parameter-set ID. Higher priorities therefore denote maintainer-curated automatic preference,
not a universal scientific quality ranking. The result retains applicability diagnostics when no
candidate can be calculated.

The current `chargefw` executable is a **water demonstration**. It builds two water conformers in
code, loads bundled parameter sets, autodetects the highest-priority applicable method and parameter
set, and prints charges. It is not yet a user-facing file/SMILES CLI.

### Calculation granularity

`Method::calculate()` operates on one molecule and zero or one conformer. The collection-level
`methods::calculate_charges()` operation applies one selected method/parameter candidate across the
prepared collection: it calculates every conformer of every molecule for geometry-dependent
methods, and one conformer-independent assignment per molecule otherwise. A molecule without a
conformer remains valid for methods that do not require geometry; geometry-dependent methods report
it as inapplicable. Future streaming and batch helpers must preserve molecule, conformer, and atom
order and report each target independently.

### Calculation selection and result contract

`calculation::calculate()` evaluates all supplied method and parameter-set candidates, then selects
the applicable candidate with the highest method priority, followed by parameter-set priority.
Equal priorities are resolved deterministically by method ID and then parameter-set ID in
lexicographic order. A successful result contains one `charges::ChargeSet`; its owned method ID and
optional parameter-set ID identify the candidate actually used. `CalculationResult::applicability`
retains the considered applicable and rejected candidates for diagnostics.

`ChargeSet` preserves calculation targets as `ChargeAssignment` entries. For geometry-dependent
methods, it contains one assignment for every conformer of every input molecule, each identified by
its molecule index and conformer index. For geometry-independent methods, it contains one
assignment per molecule and no conformer index. Each assignment's atomic-charge vector follows the
source molecule's atom order. If no candidate is applicable, `CalculationResult` contains no
`ChargeSet` and retains applicability diagnostics; calculation failures after selection are reported
as failures rather than silently treated as inapplicability.

`CalculationRequest` is intentionally a low-level native view: it contains a prepared collection,
method pointers, and spans over caller-owned parameter data. `ApplicationCalculationRequest` is the
binding-friendly owned facade: it accepts a native molecule collection, owns parameter sets, and
selects registered methods and supplied parameter sets by optional IDs. Omitted IDs use the same
deterministic automatic selection as the low-level path; explicit unavailable or inapplicable IDs
fail rather than silently falling back. Method-specific options remain a low-level advanced-native
feature until their application-facing policy is specified. All integrations must compose these
calculation paths rather than reimplement selection or scientific behavior.

### Scalable execution policy

Large-molecule acceleration is a primary product capability, not deferred work. The application
calculation contract should expose an execution policy distinct from scientific method options:
`full` performs the exact/reference calculation, while `cutoff(radius)` and `cover(radius)` perform
explicitly reported fragment approximations. Keeping execution policy outside individual `Method`
implementations allows the same spatial decomposition, fragment mapping, parallel scheduling, and
charge-reconciliation machinery to be reused by EEM/QEq-like and SQE-family methods.

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
simulation workflows possible without package-specific bindings. Biopython is a useful secondary
adapter for structural-biology object models, but coordinates and hierarchy do not guarantee a
complete chemical graph or bond orders; its connectivity, model, component, and alternate-location
policies must therefore be explicit. Gemmi is the preferred optional native C++ candidate for
PDB/mmCIF handling when a structural-biology consumer requires it.

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
broad SMILES/SDF/Mol2 workflows; optional Gemmi/Biopython adapters cover structural biology. Each
adapter translates to the same toolkit-neutral molecule and calculation contracts, so adding a
format or package cannot introduce a separate selection or scientific-policy implementation.

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
