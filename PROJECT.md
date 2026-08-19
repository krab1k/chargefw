# ChargeFW Project Guide

This document is the source of truth for implemented architecture, capabilities, compatibility state,
and product direction. See [TODO.md](TODO.md) for unfinished work, [README.md](README.md) for usage,
and [AGENTS.md](AGENTS.md) for implementation rules.

## Purpose and status

ChargeFW is a C++23, library-first framework for empirical partial atomic-charge calculation. It is a
modern successor to ChargeFW2, the engine used by Atomic Charge Calculator III (ACC III). ChargeFW is
not yet the ACC III backend; `old/` is an archived ChargeFW2 copy used for compatibility research.

The repository currently provides:

- a toolkit-neutral molecular graph and conformer model;
- prepared topology/geometry features and immutable parameter classification;
- 22 built-in methods and bundled JSON parameter sets;
- scientific applicability and execution-availability assessment;
- deterministic method/parameter selection and owned application facades;
- exact full execution and serial spatial cutoff for eight validated methods;
- native MOL/SDF/MOL2/JSON and Gemmi-backed PDB/mmCIF input;
- JSON, SDF, MOL2, and mmCIF charge output through a focused CLI.

Cover execution, parallel cutoff, broad compatibility/accuracy validation, Python bindings, packaged
distribution, and a stable failure-capable result schema remain unfinished.

## Architecture

```text
input adapters / native callers
            |
            v
core::MoleculeCollection
  atoms + bonds + zero or more conformers
            |
            v
features::PreparedMoleculeCollection
  cached topology; conformer geometry; spatial fragments
            |
            +-----------------------------+
            |                             |
            v                             v
parameters::ParameterSet          methods::MethodRegistry
            |                             |
            v                             v
ParameterClassification       MethodRequirements/options
            +-------------+---------------+
                          v
methods::find_applicable_methods()
  scientific candidates + stored classifications + full/cutoff/cover assessments
                          |
                          v
calculation::select_execution_plan()
  deterministic candidate + concrete ExecutionPolicy
                          |
                          v
calculation::calculate()
  full method execution or shared serial cutoff executor
                          |
                          v
charges::ChargeSet + application/export provenance
```

### Layer boundaries

- `core` owns only toolkit-neutral graph and conformer data.
- `features` owns derived topology/geometry state and `SpatialFragment`.
- `parameters` owns immutable matching data and classification.
- `methods` owns stateless algorithms, requirements, applicability, and the registry.
- `calculation` owns resource policy, deterministic planning, full/cutoff dispatch, and facades.
- `adapters` translate representations and preserve source identity/mapping; they do not select methods
  or duplicate scientific policy.

## Calculation contracts

### Explicit native workflow

```text
PreparedMoleculeCollection
  + ApplicabilityRequest {
      methods, parameter sets, classification options, resource policy
    }
  -> ApplicabilityResult {
      applicable candidates with classifications and execution assessments,
      rejected scientific candidates
    }
  -> caller choice or select_execution_plan()
  -> CalculationRequest {
      prepared molecules, selected ApplicableMethod, concrete ExecutionPolicy
    }
  -> CalculationResult { ChargeSet }
```

`ApplicableMethod` stores the resolved atom/bond classifications and default method options.
`CalculationRequest` does not accept classification policy and never reclassifies. The selected
candidate and prepared collection are non-owning low-level views.

`select_applicable_method()` remains available when a caller only needs priority-based scientific
candidate selection. `select_execution_plan()` additionally resolves a concrete execution mode.

### Owned application workflow

`ApplicationCalculationRequest` owns molecules and parameter sets and accepts optional method and
parameter-set IDs, classification options, an `ExecutionSelection`, and a `ResourcePolicy`.

- `assess()` prepares molecules, finds applicable candidates, and selects a concrete plan without
  calculating.
- `calculate()` composes assessment and execution.
- Omitted IDs use deterministic ranking: method priority, parameter-set priority, method ID, then
  parameter-set ID.
- Explicit unavailable/inapplicable IDs or unsupported explicit execution fail; there is no fallback.
- Application-facing method-option selection is not implemented; applicability currently uses each
  method's validated defaults.

`ChargeSet` owns the selected method ID and optional parameter-set ID. Geometry-dependent methods
produce one source-ordered assignment per molecule conformer; geometry-independent methods produce one
assignment per molecule.

## Execution policy and cutoff

`ExecutionSelection` represents caller preference (`auto`, `full`, `cutoff`, or `cover`). Results use
a concrete `ExecutionPolicy`; `auto` never survives as the effective mode.

- Full execution accepts no radius or charge correction.
- Reduced radii must be finite and at least 8 Å.
- Explicit cutoff/cover requires a radius.
- Automatic reduced execution uses 12 Å unless the caller supplies a valid radius.
- Reduced execution defaults to uniform final charge correction; explicit cutoff may select `none`.
- The shared full atom threshold defaults to 20,000; `nullopt`/`unlimited` disables only this safeguard.
- A method is considered expensive when declared cubic in time or quadratic in memory. Above a finite
  threshold, automatic planning avoids expensive full execution; explicit full remains allowed and
  returns a resource warning.
- Threshold assessment is per molecule, while one execution mode applies to the collection.

Scientific applicability is independent of the resource warning. Missing coordinates, topology,
formal charges, element properties, required parameter coverage, or other method prerequisites remain
hard failures.

### Implemented cutoff methods

Serial cutoff is declared for:

```text
abeem, eem, eqeq, eqeqc, qeq, sqe, sqeq0, sqeqp
```

The shared executor builds one source-ordered induced spatial fragment per source atom, projects the
whole-molecule classification, invokes the selected method through ordinary `CalculationInput`, and
keeps the explicitly mapped center charge. It then applies the selected final correction. Neighbor
search is currently linear and execution is serial.

EEM/QEq-like methods and ABEEM allocate each fragment a formal-charge target proportional to its atom
count and restore the source formal-charge total. SQE uses zero fragment and final targets. SQE+q0
projects formal initial charges; SQE+qp projects parameterized `q0` corrected to fragment formal
charge. SQE-family cutoff is a new approximation design, not ChargeFW2 parity.

Whole-molecule-radius tests verify full/cutoff agreement for all eight methods on small fixtures,
including neutral/charged and zero-width SQE-family cases. They do not establish general cutoff error
envelopes. Cover is represented in policy and applicability but no method declares it supported and
its executor is not implemented.

### Preliminary large-structure observations

Manual development runs, not supported benchmarks, show bounded fragment memory and convergence toward
full execution:

- `10aw.cif` (10,479 atoms): QEq cutoff MAE decreased from 0.004288 e at 8 Å to 0.001751 e at 12 Å;
  cutoff peak RSS was about 30 MB versus 1.76 GB for full.
- At 10 Å on `10aw.cif`, EEM, QEq, and EQeq produced finite source-ordered neutral assignments with
  MAE of 0.001689 e, 0.002521 e, and 0.000379 e, respectively.
- At 10 Å on `1ek9.cif` polymers (9,798 atoms), all eight supported methods completed in full and
  cutoff modes. SQE+q0 had materially larger deviation (0.053125 e MAE) than the other tested methods,
  so it needs charged/disconnected and multi-radius validation before an accuracy claim.
- SFKEEM is intentionally not cutoff-capable: local neutral-fragment target charges lose its global
  chemical-potential constraint. It needs a dedicated truncated-kernel solver or a separately
  validated fragment policy.

Hardware, build, parameter data, charge state, and structure affect these observations. They are
regression-investigation references only, not automatic-policy evidence or compatibility tolerances.

See [CUTOFF_IMPLEMENTATION_PLAN.md](CUTOFF_IMPLEMENTATION_PLAN.md) for the retained design decisions
and validation still required.

## Methods and parameters

The registry contains 22 method IDs:

```text
abeem, charge2, delre, denr, dummy, eem, eqeq, eqeqc, formal, gdac,
kcm, mgc, mpeoe, peoe, qeq, sfkeem, smpqeq, sqe, sqeq0, sqeqp, tsef, veem
```

Bundled parameter sets live in `data/parameters/` and install under
`share/chargefw/parameters`. All nine archived SQE-family sets are migrated. Numerical parity with
ChargeFW2 remains incomplete and release-blocking for production adoption.

## Molecular I/O and CLI

The `chargefw` executable has `calculate`, `inspect`, `applicability`, `methods`, and `parameters`
subcommands. It selects input by extension:

| Input | Reader | Calculation output |
| --- | --- | --- |
| `.mol`, `.sdf` | Native MOL/SDF | JSON, SDF, MOL2, mmCIF |
| `.mol2` | Native MOL2 | JSON, SDF, MOL2, mmCIF |
| `.json` | ChargeFW schema 1.0 | JSON, SDF, MOL2, mmCIF |
| `.pdb` | Gemmi PDB | JSON, mmCIF |
| `.cif`, `.mmcif` | Gemmi mmCIF | JSON, mmCIF |

The CLI currently reads the complete collection before calculation and rejects the input on the first
malformed record. It is not yet a bounded-memory batch runner. JSON input with multiple conformers is
rejected by molecular output because SDF/MOL2 export currently accepts one assignment per molecule.

Structural input supports record selection (`all`, `polymers-and-ligands`, `polymers`) and connectivity
(`none`, `explicit`, `templates`, `hybrid`). Library structural readers default to `none`; the CLI
deliberately defaults to `hybrid`. Alternate locations select blank, then `A`, then the first
occurrence. Compact templates cover standard amino acids, standard RNA/DNA nucleotides, water, and
sequential peptide/nucleotide links. There is no distance-based bond perception or full CCD provider.

Adapters preserve selected source atom order, formal charges, conformer identity, and record identity.
MOL/SDF supports a deliberately narrow V2000/V3000 subset; MOL2 accepts standard element-prefixed atom
types and numeric bonds. Aromatic bonds are imported as single bonds. Partial charges in MOL2 input are
ignored rather than treated as formal charges.

SDF and MOL2 same-format writers preserve source content while adding/replacing ChargeFW-owned charge
fields. The Gemmi writer semantically preserves mmCIF categories (presentation may normalize), converts
PDB through Gemmi, and generates local `UNL` blocks for nonstructural input.

### JSON result state

Schema `1.0` records source identity, source-to-native atom/conformer mapping, source-ordered
assignments, totals, and diagnostics. Invocation-level `calculation_provenance` records requested
method/parameter IDs, permissive typing, resource threshold, execution/radius/correction, structural
input policy, and effective method/parameter/execution plus warnings. Charges are rounded to at most
four decimal places in JSON only.

Current limitations:

- unsuccessful automatic selection emits a JSON error document to stdout rather than the requested
  output directory;
- import and calculation exceptions are not represented as owned record-scoped result entries;
- exit statuses do not distinguish invalid requests from calculation failures;
- SDF, MOL2, and mmCIF do not yet carry the complete JSON provenance;
- explicit-full threshold warnings are retained in the result but are not emitted before solver
  allocation.

## Build, installation, and distribution

CMake builds one `chargefw_core` library and an optional `chargefw` CLI. Public dependencies are
nlohmann/json and Gemmi; Eigen is private and CLI11 is used by the application. CMake first searches
for compatible packages and otherwise fetches CLI11 2.6.2, nlohmann/json 3.12.0, Eigen 5.0.1, and Gemmi
0.7.4.

Installation currently provides the library, public headers, CLI, generated config header, and bundled
parameter JSON. Exported CMake package targets are not implemented. There are no repository CI
workflows, Python bindings/wheels, Conda recipe, or container image yet.

## Integration direction

The intended front ends converge on the same owned application facade:

```text
native files / C++ / future NumPy and RDKit converters
                         |
                         v
toolkit-neutral molecule data -> native assessment/calculation -> charges + mapping + provenance
```

The first Python API should use NumPy-style arrays and nanobind. RDKit integration should initially be
a pure-Python converter with lazy optional import, preserving atom/conformer identity and performing no
implicit sanitization, hydrogen changes, protonation, embedding, or optimization. Packaging work must
not make RDKit a dependency of the core library or base wheel.

## Scientific and compatibility principles

- Preserve source atom order and identify every molecule/conformer assignment.
- Record method, parameter set, options when exposed, execution approximation, correction, and warnings.
- Keep exact full calculations distinct from cutoff/cover approximations.
- Make automatic policy explainable and explicitly overridable.
- Compare ChargeFW2 with identical graphs, charges, coordinates, options, and parameter data.
- Record per-method tolerances and intentional deviations in tests or maintained compatibility data.
- Do not claim cutoff accuracy from whole-fragment equality or one-off large-structure measurements.

The prioritized unfinished work and acceptance criteria are maintained only in [TODO.md](TODO.md).
