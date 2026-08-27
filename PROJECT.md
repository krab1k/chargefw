# ChargeFW Project Guide

This document is the source of truth for implemented architecture, capabilities, compatibility state,
and product direction. See [TODO.md](TODO.md) for unfinished work, [PYTHON.md](PYTHON.md) for the Python
API and delivery plan, [README.md](README.md) for usage, and [AGENTS.md](AGENTS.md) for implementation
rules.

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
- focused planning and execution contracts for no-plan, explicit-policy, ownership, cardinality,
  source-order, and reduced failure behavior;
- exact full execution, oneTBB-parallelized KD-tree spatial cutoff and spatial cover for nine
  reduced-capable methods;
- native MOL/SDF/MOL2/JSON and Gemmi-backed PDB/mmCIF input;
- JSON, SDF, MOL2, and mmCIF charge output through a focused CLI.

Independent broad-accuracy and reduced-mode validation, relocatable native packaging, Python bindings,
and packaged distribution remain unfinished.

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
  mode dispatch -> full/cutoff/cover executor -> shared target execution
                          |
                          v
charges::ChargeSet + application/export provenance
```

### Layer boundaries

- `core` owns only toolkit-neutral graph and conformer data.
- `features` owns derived topology/geometry state and `SpatialFragment`.
- `parameters` owns immutable matching data and classification.
- `methods` owns stateless algorithms, requirements, applicability, and the registry.
- `calculation` owns resource policy, deterministic planning, shared full/cutoff/cover dispatch, and
  facades.
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
candidate selection. The owned assessment expands applicable candidates into concrete execution plans.

### Owned application workflow

`AssessmentRequest` owns molecules and parameter sets and accepts optional method and
parameter-set IDs, classification options, an `ExecutionSelection`, and a `ResourcePolicy`.
Parameter-set IDs must be unique across an assessment request, including sets for different methods;
duplicate IDs are rejected before parameter filtering and applicability assessment.

- `assess(request)` copies an lvalue request's molecule and parameter inputs; `assess(std::move(request))`
  snapshots its selection configuration, then transfers its execution inputs. Callers must not inspect a
  request after it is consumed. Both forms prepare molecules, find applicable candidates, and produce
  priority-ordered concrete plans plus rejected policy alternatives without calculating. Their result
  keeps prepared features, classifications, methods, and parameter lifetimes private.
- `calculate(assessment, plan, max_threads, observer)` executes an exact target-bound plan without
  repeating preparation or classification. Assessments and plans remain reusable, including for
  independent concurrent executions. The default-plan convenience overload returns a no-plan result when
  no plan exists. `ExecutionResult::effective` records detached method, parameter-set, validated-option,
  concrete-policy, and warning provenance. Thread limits and observation are execution-only inputs.
- Rejected application-facing candidates retain the method identity from the selected applicability
  method list, including explicitly requested methods, rather than using registry ordering.
- Omitted IDs use deterministic ranking: method priority, parameter-set priority, method ID, then
  parameter-set ID.
- Explicit unavailable/inapplicable IDs or unsupported explicit execution fail; there is no fallback.
- Application-facing method options are supplied as method-scoped overrides. Explicit method selection
  rejects options for other methods; automatic selection accepts namespaced overrides and records the
  selected method's complete validated options in result provenance.

`ChargeSet` owns the selected method ID and optional parameter-set ID. Geometry-dependent methods
produce one source-ordered assignment per molecule conformer; geometry-independent methods produce one
assignment per molecule.

## Execution policy and cutoff

`ExecutionSelection` represents caller preference (`auto`, `full`, `cutoff`, or `cover`). Results use
a concrete `ExecutionPolicy`; `auto` never survives as the effective mode.
Domain-owned `*_from_string` and `to_string` functions define the canonical execution, method-report,
and structural-input vocabulary shared by the CLI and Python binding; frontends do not duplicate enum
switches or invent alternate spellings.

- Full execution accepts no radius or charge correction.
- Reduced radii must be finite and at least 8 Å.
- Explicit cutoff/cover requires a radius.
- Automatic reduced execution uses 12 Å unless the caller supplies a valid radius.
- Reduced execution defaults to uniform final charge correction; explicit cutoff may select `none`.
- The automatic cutoff threshold defaults to 20,000 atoms. For expensive full methods, exceeding it
  promotes automatic planning to cutoff; `nullopt`/`unlimited` disables this safeguard.
- The automatic cover threshold defaults to 80,000 atoms. Exceeding it promotes automatic cutoff
  planning to cover; `nullopt`/`unlimited` keeps cutoff eligible without a size limit. A finite cover
  threshold requires a finite cutoff threshold and must not be smaller.
- A method is considered expensive when declared cubic in time or quadratic in memory. Above a finite
  threshold, automatic planning avoids expensive full execution; explicit full remains allowed and
  returns a resource warning.
- Threshold assessment is per molecule, while one execution mode applies to the collection.

Scientific applicability is independent of the resource warning. Missing coordinates, topology,
formal charges, element properties, required parameter coverage, or other method prerequisites remain
hard failures.

### Calculation observability

Calculations report progress and support cooperative cancellation through a per-request
`CalculationObserver` (`include/chargefw/calculation/observer.h`). The observer is non-owning:
callers retain the object and pass a reference through `CalculationRequest` or `calculate()`. Calls
that do not need progress use the stateless `default_calculation_observer()`; execution always uses
the same observed loop.

Events are structured as `CalculationProgress` snapshots carrying the execution mode, method ID,
target index/count, completed-fragment count/total, molecule/conformer identity, and elapsed seconds.
Target indices are zero-based; completed-fragment counts start at one and are not fragment identifiers.
Non-owning views in the snapshot are valid only during the callback.

Three event tiers exist:

- **Computation phases** — `computation_started/finished`. The low-level
  `calculate(const CalculationRequest&)` owns this tier, so both public calculation overloads emit
  each start with exactly one terminal `computation_finished` event on success, cancellation, and
  non-cancellation calculation failure. The low-level boundary propagates cancellation and failures;
  the owned facade converts cancellation to its result. Both events carry the effective mode and
  method ID. Assessment is not observed.
- **Per-target events** — `target_started/finished` for each molecule+conformer pair. Emitted by the
  executors with correct per-target `molecule_index` and `conformer_index`.
- **Fragment progress** — `fragment_progress` reports aggregate completed/total counts for cutoff
  (center-atom fragments) and cover (pivot fragments). The first completed fragment is reported
  immediately; later snapshots are throttled to one per 200 ms. Counts strictly increase per target
  among emitted snapshots but may skip completed fragments. A non-empty successful target has exactly
  one terminal `count/count` snapshot; empty targets have none. Parallel callbacks for different
  targets may run concurrently, so global callback arrival order does not establish completion order.

Per-plan warnings, including explicit-full execution above the resource threshold, are available through
`ExecutionPlan::warnings()` before execution begins. Callers may report them before calling
`calculate(assessment, plan, ...)`; they are also retained in detached effective result provenance.

Cooperative cancellation: `cancelled()` is polled at each target and fragment iteration boundary,
including immediately after `target_started`. When it returns true, `CalculationCancelled` is thrown
and propagated through oneTBB to the application facade, which returns
`ExecutionResult{.status = ExecutionStatus::cancelled}` with no charges. Low-level `calculate()`
propagates the exception. The observer contract requires its callbacks and `cancelled` to be
thread-safe and must not throw; the observer is purely observational and must not mutate method
options, parameters, execution policy, geometry, or selection. This terminal-event guarantee lets
terminal observers restore progress output on all started-computation paths.

The owned facade represents `success`, `invalid_input_or_request`, `no_executable_plan`,
`numerical_failure`, and `cancelled` outcomes. The CLI writes a JSON result for all outcomes reached
after import, with document and source-record status/diagnostics; non-successful outcomes never write
charge assignments or molecular output. CLI exits are `0` for success, `1` for an unexpected internal
failure, `2` for invalid input/request, `3` for no executable plan, `4` for numerical failure, and `5`
for cancellation.

Target calculation failures retain method, one-based molecule, molecule-name, and conformer context.
The owned facade preserves invalid request exceptions for the caller/CLI input boundary, converts other
method or solver failures into a `numerical_failure` result with that context, and preserves cancellation.
The singular EEM regression exercises this path without treating it as a general conditioning benchmark.

The CLI progress observer serializes standard-error updates and clears the full rendered line before
every target/fragment update and terminal newline, preventing stale suffixes when messages shrink or
their display tier changes.

Facade and observer regression coverage includes rejected-report identity, explicit no-plan and
unsupported-policy boundaries, empty/tiny input cardinality, full and reduced solver failures,
fragment-progress cancellation in serial and parallel cutoff/cover execution, ownership after
lvalue/rvalue assessment, and source-ordered multi-molecule/multi-conformer results.

### Implemented reduced methods

Cutoff and cover are declared for:

```text
abeem, eem, eqeq, eqeqc, qeq, sfkeem, sqe, sqeq0, sqeqp
```

The shared executor builds one source-ordered induced spatial fragment per source atom, projects the
whole-molecule classification, invokes the selected method through ordinary `CalculationInput`, and
keeps the explicitly mapped center charge. It then applies the selected final correction. Neighbor
search uses one `SpatialFragmentBuilder` and nanoflann KD-tree per source conformer; execution is
parallel over independent molecule/conformer targets for full execution. Cutoff and cover process
targets serially and parallelize the active target's fragment centers or pivots with the caller's full
worker budget. Nested scheduling is avoided and result materialization remains source-ordered.

EEM/QEq-like methods and ABEEM allocate each fragment a formal-charge target proportional to its atom
count and restore the source formal-charge total. SQE uses zero fragment and final targets. SQE+q0
projects formal initial charges; SQE+qp projects parameterized `q0` corrected to fragment formal
charge. SQE-family cutoff is a new approximation design, not ChargeFW2 parity.

Cover constructs source-order pivots and ownership serially, then builds and solves each independent
radius fragment in parallel. It retains charges for every previously unassigned atom within 3 Å of the
pivot. Overlapping retained interiors use first-pivot ownership; halo-only atoms provide context.
Cover shares classification projection, fragment target charges, validation, and final correction with
cutoff. Its pivot fragments use the reduced execution policy above; result materialization remains
source-ordered.

All execution modes validate selected-candidate parameter/classification invariants and reject
coordinate-dependent candidates for molecules without conformers before calculation.
Method-neutral regression tests cover serial pivot ownership, retained-interior coverage and first-pivot
overlap reconciliation, fragment/whole target charges and correction, and bit-identical serial versus
parallel source-ordered output.

Whole-molecule-radius tests verify full/cutoff/cover agreement for all eight methods on small fixtures,
including neutral/charged and zero-width SQE-family cases. They do not establish general reduced-mode
error envelopes.

### Preliminary large-structure observations

Manual development runs and local query measurements, not supported benchmarks, show bounded fragment
memory and convergence toward full execution:

- `10aw.cif` (10,479 atoms): QEq cutoff MAE decreased from 0.004288 e at 8 Å to 0.001751 e at 12 Å;
  cutoff peak RSS was about 30 MB versus 1.76 GB for full.
- At 10 Å on `10aw.cif`, EEM, QEq, and EQeq produced finite source-ordered neutral assignments with
  MAE of 0.001689 e, 0.002521 e, and 0.000379 e, respectively.
- At 10 Å on `1ek9.cif` polymers (9,798 atoms), all eight supported methods completed in full and
  cutoff modes. SQE+q0 had materially larger deviation (0.053125 e MAE) than the other tested methods,
  so it needs charged/disconnected and multi-radius validation before an accuracy claim.
- A GCC release query comparison over all centers found the KD-tree radius search 1.57–3.99× faster at
  12 Å and 3.96–9.70× faster at 8 Å than the former linear scan on `10aw.cif` (10,479 atoms),
  `1ek9.cif` (11,426 atoms), and `8yax.cif` (40,738 atoms). Per-conformer tree construction took
  0.8–3.8 ms; source-order result materialization limits the advantage as radius and neighbor count
  grow.

Hardware, build, parameter data, charge state, and structure affect these observations. They are
regression-investigation references only, not automatic-policy evidence or compatibility tolerances.

Automatic planning examines candidates in deterministic method/parameter priority order and selects,
per candidate, non-discouraged full execution, cutoff below the cover threshold, then cover. It never
silently forces discouraged full or cutoff execution. Unfinished reduced-execution work and its
validation gates are tracked in [TODO.md](TODO.md).

## Methods and parameters

The registry contains 22 method IDs:

```text
abeem, charge2, delre, denr, dummy, eem, eqeq, eqeqc, formal, gdac,
kcm, mgc, mpeoe, peoe, qeq, sfkeem, smpqeq, sqe, sqeq0, sqeqp, tsef, veem
```

Bundled parameter sets live in `data/parameters/` and install under
`share/chargefw/parameters`. All nine archived SQE-family sets are migrated. Scientific comparison
against ChargeFW2 and ACC III remains incomplete and release-blocking for production adoption; parity
with the archived implementation is not assumed where the publication, inputs, or current design
justify a difference.

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
malformed record. It is not yet a bounded-memory batch runner. `--conformers first|all` (default `all`)
selects the first or all conformers/models for every reader and is ignored for molecules without
conformers. JSON input with multiple conformers is rejected by SDF/MOL2 output because those formats
currently accept one assignment per molecule; `--conformers first` can select a compatible output.
Structural output may retain uncalculated source models when `first` is selected, while charge rows
are emitted only for calculated conformers and use the corresponding source atom IDs. The CLI partitions
imported execution molecules from retained source/export context before assessment. It snapshots requested
provenance and execution thread policy before consuming the assessment request; applicability-only
execution transfers imported molecules because it has no export phase.

Generated SDF and MOL2 output uses source conformer zero for a geometry-independent assignment. It
therefore requires source coordinates but does not require the selected method to depend on geometry.

Structural input supports record selection (`all`, `polymers-and-ligands`, `polymers`) and connectivity
(`none`, `explicit`, `templates`, `hybrid`). Library structural readers default to `none`; the CLI
deliberately defaults to `hybrid`. Alternate locations select blank, then `A`, then the first
occurrence. Compact templates cover standard amino acids, standard RNA/DNA nucleotides, water, and
sequential peptide/nucleotide links. Sequential links require the same chain and consecutive residue
sequence numbers, with insertion-code continuity allowed when the first inserted residue keeps the
parent sequence number. There is no distance-based bond perception or full CCD provider. Structural
readers reject empty models and conformer sequences whose selected atom identity, element, formal
charge, or name differ across models.

Adapters preserve selected source atom order, formal charges, selected conformer identity, and record
identity. For mmCIF, the original document and imported block mapping are retained for
preservation-oriented output: unselected alternate-location rows remain in the source structure, while
charge rows reference only the selected `_atom_site.id` values. PDB input is converted to normalized
mmCIF containing only selected atoms.

Calculation JSON includes execution metrics for parsing, applicability, computation, non-JSON output
writing, wall-clock timestamps, and peak resident memory. Applicability and computation timings use
monotonic clocks; timestamps are UTC. Feature preparation occurs during assessment and is retained for
execution, so applicability timing includes preparation and computation timing covers execution only.
MOL/SDF supports a deliberately narrow V2000/V3000 subset; MOL2 accepts standard element-prefixed atom
types and numeric bonds. Aromatic bonds are imported as single bonds. Partial charges in MOL2 input are
ignored rather than treated as formal charges.

SDF and MOL2 same-format writers preserve source content while adding/replacing ChargeFW-owned charge
fields. Each SDF charge vector is identified by property suffix and source atom order and has a compact
`CHARGEFW_CHARGE_METADATA_N` field: `type`, `method`, `parameter_set` (`.` when inapplicable),
`software_name`, and `software_version`. mmCIF uses the same compact assignment-level attribution in
the published `mmcif_charges_v11.dic`, with dictionary assignment IDs and source atom IDs. Complete
invocation provenance remains in the companion JSON result. MOL2 has no portable provenance extension,
so its output is intentionally limited to the atom partial-charge column and preserved source content.
The Gemmi writer semantically preserves mmCIF categories (presentation may normalize), converts PDB
through Gemmi, and generates local `UNL` blocks for nonstructural input. Its universal partial-charge
dictionary attribution is not intended to duplicate complete ChargeFW invocation provenance.

### JSON result state

Schema `1.0` records source identity, source-ordered assignments, totals, and diagnostics. Input parsing
remains fail-fast: malformed input rejects the document rather than manufacturing partial failed records.
Once import and request construction succeed, invalid assessment, no-plan, numerical-failure, and
cancellation outcomes use the same owned result document and status vocabulary. Successful adapter
warnings, such as ignored MOL2 partial charges, remain attached to their source record.

Diagnostics have a constrained severity, stable code, explanatory message, and optional zero-based
molecule, atom, bond, and conformer indices or one-based source line. Human-readable messages use
one-based molecule, atom, bond, and conformer numbering and include available chemical context such as
element, source atom name, formal charge, bond participants, and bond order. No-plan results retain the
actual rejected-candidate reasons instead of only reporting a generic selection failure.

Geometry-dependent methods write one assignment per selected conformer; geometry-independent methods
write one assignment per molecule. Charges are rounded to at most four decimal places in JSON only, and
each `total_charge` is the sum of those serialized values. Internal calculations retain native precision.
The result schema remains `1.0` while the project is pre-release and may be refined in place. After the
first stable release, incompatible changes require a major schema-version change; additive optional
fields require a minor change and must preserve the meaning of existing fields.
`generator.version` identifies the ChargeFW release, including its method implementation,
classification behavior, and bundled publication-derived parameter data. Effective parameter-set IDs
identify the selected bundled set within that release. Reproducible results require matching software
version, result-schema version, effective method and parameter-set IDs, input/coordinates, options, and
execution policy. Custom caller-supplied parameter sets are outside this bundled-data identity guarantee.
Invocation-level `calculation_provenance` records requested conformer selection, method/parameter IDs,
permissive typing, cutoff/cover resource thresholds, execution/radius/correction, structural input
policy, and effective method/parameter/execution plus warnings.

Current limitations:

- failures before source import or request construction completes cannot produce a source-record result
  document;
- MOL2 and mmCIF intentionally do not carry complete JSON provenance: MOL2 carries charge values
  without portable provenance metadata, while mmCIF uses the universal partial-charge dictionary.

## Build, installation, and distribution

CMake builds one `chargefw_core` library and an optional `chargefw` CLI. Public dependencies are
nlohmann/json and Gemmi; Eigen, nanoflann, and oneTBB are private and CLI11 is used by the application.
CMake fetches pinned CLI11 2.7.2, nlohmann/json 3.12.0, Eigen 5.0.1, nanoflann 1.12.1, oneTBB 2023.1.0,
and Gemmi 0.7.4 by default. `CHARGEFW_USE_SYSTEM_DEPENDENCIES=ON` opts into find-package-first behavior.
The CLI and test suite default to off when ChargeFW is configured as a subproject.

Installation provides the library, public headers, CLI, generated config header, bundled parameter JSON,
and an exported `chargefw::core` CMake target. Default installations include nlohmann/json and Gemmi
development packages and the private oneTBB runtime; private build-only dependencies and test tooling are
excluded. There are no repository CI workflows, Python bindings/wheels, Conda recipe, or container image yet.

Default parameter discovery resolves installed JSON relative to the loaded ChargeFW library, so an
installed prefix remains usable after it is moved. Build-tree CLI execution is unsupported; CLI tests
install to a temporary prefix before running. A downstream CMake smoke test configures, builds, and runs
against the moved prefix using `find_package(chargefw CONFIG REQUIRED)`. A separate smoke test moves and
runs a custom `lib/chargefw` plus `resources` installation layout.

## Integration direction

The intended front ends converge on the same owned application facade:

```text
native files / C++ / future NumPy and RDKit converters
                         |
                         v
toolkit-neutral molecule data -> native assessment/calculation -> charges + mapping + provenance
```

The first Python API should use NumPy-style arrays and nanobind. Its detailed public contract, required
Gemmi integration, optional toolkit boundaries, packaging design, and phased delivery plan are maintained
in [PYTHON.md](PYTHON.md). Packaging work must not make RDKit or Biopython a dependency of the core
library or base wheel.

## Scientific and compatibility principles

### ChargeFW2 and publication audit

The bounded audit covered the 20 scientific methods shared with ChargeFW2, their bundled parameter sets
where applicable, a canonical CHNO control, ionic and geometry-independent controls, and a compact
representative panel. Current and legacy calculations have exact or near-machine numerical parity under
identical finite-molecule inputs and algorithm-parity options. Methods whose cited neutral-only
formulation cannot satisfy ionic requests reject those inputs.

Publication review supports the implemented finite, non-periodic variants and records intentional scope
differences. In particular, periodic/Ewald branches, EQeq charge centers, and original-QEq Slater,
hydrogen-iteration, and charge-bound branches are unsupported. The audit is compatibility and
conformance evidence, not a general QM-accuracy or reduced-mode-accuracy claim.

The committed representative execution panel covers water, ethanol, benzene, acetate, methylammonium,
glycine zwitterion, dimethyl sulfoxide, chlorobenzene, methyl phosphate, and heme. It provides stable,
inspectable behavior controls for calculation and full/cutoff/cover tests; it is not a substitute for a
method-specific reduced-mode accuracy study.

- Preserve source atom order and identify every molecule/conformer assignment.
- Record method, parameter set, options when exposed, execution approximation, correction, and warnings.
- Keep exact full calculations distinct from cutoff/cover approximations.
- Make automatic policy explainable and explicitly overridable.
- Compare ChargeFW2 with identical graphs, charges, coordinates, options, and parameter data, while
  treating its output as evidence for investigation rather than an automatic oracle.
- Investigate and classify per-method differences as a current defect, a legacy defect, an intentional
  change, or unresolved; record tolerances and rationale in tests or maintained comparison data.
- Do not claim cutoff accuracy from whole-fragment equality or one-off large-structure measurements.

The prioritized unfinished work and acceptance criteria are maintained only in [TODO.md](TODO.md).
