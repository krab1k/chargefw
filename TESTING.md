# ChargeFW Testing Strategy

This document defines how the repository's tests should be organized and improved. It is a working
plan for implementation-sanity, contract, and integration testing. It does not replace the scientific
and compatibility evidence tracked in [TODO.md](TODO.md) or the architecture in [PROJECT.md](PROJECT.md).

## Goal

Make the test suite sensible as a whole: each test should have a clear contract, the appropriate
layer, a focused fixture, and a non-duplicated reason to exist. Adding a feature should extend an
existing test category where possible rather than creating another broad catch-all test.

The suite must establish that the implementation obeys ChargeFW's declared behavior. It is not, by
itself, a claim that a method is scientifically accurate or compatible with ChargeFW2/ACC III.
Scientific comparisons, reduced-mode error envelopes, and published-reference validation remain
separate evidence work.

## Testing principles

### 1. Implementation sanity is valuable and is not scientific validation

A small chemically intelligible fixture can expose ordinary implementation errors without being a
scientific benchmark. Neutral water is the canonical example:

- for every charge-calculating built-in method except `formal` and `dummy`, oxygen must be negative;
- both hydrogens must be positive and equal within a method-appropriate numerical tolerance;
- the assignment must have the expected atom count and conserve the molecule's formal charge;
- geometry-dependent methods must produce one assignment for each conformer and react to a changed
  conformer where their formula depends on geometry;
- topology-only methods must produce one conformer-independent assignment and must not change when
  coordinates are added or changed.

These are implementation-sanity contracts. They detect swapped atom mappings, lost constraints,
incorrect parameter classification, wrong conformer dispatch, accidental index dependence, and many
sign or solver regressions. They do not establish a method's publication-level accuracy.

`formal` and `dummy` have their own exact contracts: formal charges are copied exactly, and dummy
charges are zero exactly.

### 2. Test the public contract at the lowest useful layer

A failure should point to one subsystem whenever practical:

- value and validation rules belong in `core`, `charges`, and policy tests;
- derived graph/geometry/mapping behavior belongs in `features`;
- matching and immutable parameter access belong in `parameters`;
- method prerequisites and algorithms belong in `methods`;
- candidate ranking, policy, target fan-out, and mode dispatch belong in `calculation`;
- representation mapping and preservation belong in `adapters`;
- process exit status, files, serialization, and terminal rendering belong in CLI integration tests.

Do not use a CLI test to prove a classifier rule, or a method test to prove execution-plan ranking.

### 3. Test contracts and invariants before reference values

Every method and execution mode needs structural checks: finite values, source order, cardinality,
identity, charge conservation where applicable, and deterministic behavior. A small number of exact
or tolerance-based numerical regression values should supplement these checks for method-specific
logic. Avoid treating a broad full-versus-reduced comparison as proof of accuracy.

### 4. Keep negative tests structured

For invalid input, inapplicability, unavailable execution, solver failure, and cancellation, assert
observable structured state and context:

- issue kind and associated method/parameter/molecule/conformer where exposed;
- whether a plan exists;
- whether charges are absent;
- result or exception category at the correct boundary;
- terminal observer behavior where calculation started.

Only assert exception wording when text itself is a documented user-facing contract. Prefer typed
issues and fields over substrings.

### 5. Checks must be active in every configured build

The approved C++ test framework is [Snitch](https://github.com/cschreib/snitch). Use named
`TEST_CASE`s and Snitch assertions such as `CHECK`, `REQUIRE`, `CHECK_THROWS_AS`, and `CHECK_THAT`;
they remain active under `NDEBUG` and report the failed contract. Do not add new standard `assert`
checks or extend `tests/cpp/support/test_assertions.h`; migrate its remaining callers progressively.

Snitch is a test-only dependency. CMake must locate it first and otherwise obtain it with
`FetchContent` only when `CHARGEFW_BUILD_TESTS` is enabled. It must not be fetched, linked by library
targets, installed, or required for non-test builds.

## Current suite: useful coverage and issues

The repository has layer-based CTest executables under `tests/cpp/`, including core, features,
parameters, methods, calculation, charges, adapters, and two CLI scripts. This is a good starting
layout. The main problems are organization and uneven depth.

### Strengths

- Core value types, molecule validation, topology, conformer geometry, spatial fragments, charge
  containers, and execution-policy value validation have focused tests.
- The calculation layer covers full/cutoff/cover dispatch, source-ordered assignments, correction,
  cancellation, observer terminal events, and several failure paths.
- Native and Gemmi adapters have meaningful reader/writer, preservation, model, and bond-strategy
  coverage.
- All built-in methods are registered and most have at least a smoke or numerical test.
- Reduced execution has whole-molecule-radius equality checks for all eight declared reduced-capable
  methods. This is a useful executor regression, not a general approximation validation.

### Structural issues

- Planning, facade behavior, target ordering, reduced execution, and observer lifecycle are partly
  duplicated across calculation tests.
- Method-specific tests still need rationalization: retain only distinguishing algorithms and stable
  numerical regressions now that generic workflow and water-sanity coverage are manifest-driven.
- Parameter classification and bundled parameter-data integrity are tested much less thoroughly than
  method execution.
- CLI tests primarily exercise successful output. Failure states, result-schema behavior, progress
  erasure, record ordering, and non-success exit distinctions remain incomplete.

### Current test inventory

Every CTest entry has one primary taxonomy layer. A target may incidentally exercise a lower-level
type, but that is not its reason to exist.

| Taxonomy | Current CTest entries |
| --- | --- |
| A — value and model | `test_atom`, `test_bond`, `test_conformer`, `test_molecule`, `test_molecule_collection`, `test_periodic_table`, `test_atomic_charges`, `test_charge_collection`, `test_execution_policy` |
| B — features and parameters | `test_topology_features`, `test_conformer_features`, `test_spatial_fragment`, `test_parameter_classifier`, `test_parameter_set_io` |
| C — built-in and method conformance | `test_builtin_methods`, `test_method_options`, `test_method_registry`, `test_method_prerequisites`, `test_method_applicability`, `test_method_calculation`, `test_calculation_input` |
| D — method-specific algorithms | `test_peoe`, `test_mpeoe`, `test_veem`, `test_gdac`, `test_tsef`, `test_delre`, `test_denr`, `test_kcm`, `test_qeq`, `test_sqe`, `test_sqeq0`, `test_sqeqp`, `test_eem`, `test_smpqeq`, `test_sfkeem`, `test_eqeq`, `test_eqeqc`, `test_abeem`, `test_mgc`, `test_charge2` |
| E — calculation planning and execution | `test_planning`, `test_calculation`, `test_observer`, `test_calculation_targets`, `test_cutoff_execution`, `test_cover_execution` |
| F — adapters and CLI | `test_mol`, `test_json`, `test_json_output`, `test_mol2_output`, `test_sdf_output`, `test_pdb`, `test_mmcif`, `test_mmcif_output`, `test_chargefw_cli_outputs`, `test_chargefw_cli_structural_outputs` |

### Known duplicate or misplaced coverage

These should be consolidated without reducing the contracts checked:

- Mixed multi-molecule/multi-conformer source ordering is asserted in full execution by
  `tests/cpp/calculation/test_calculation.cpp` and in cutoff/cover execution by
  `tests/cpp/calculation/test_cutoff_execution.cpp`. Observer tests keep only event-identity
  assertions.
- Full/cutoff/cover whole-radius equality is implemented in `test_cutoff_execution.cpp` even though it
  covers both cutoff and cover. Place it in a neutral reduced-execution contract test.

## Desired test taxonomy

### A. Value and model contracts

Small tests for public value types and validation:

- atoms, bonds, conformers, molecules, collections, periodic table;
- execution selection/policy/resource-policy validation and string conversion;
- atomic charges, assignments, charge sets, and collection validation.

These tests should be exhaustive for boundary and invalid values and should not depend on methods or
file formats.

### B. Prepared-feature and parameter contracts

Tests for derived state and immutable matching:

- topology adjacency, bond lookup, and source order;
- conformer finite/coincident-coordinate detection;
- spatial fragment membership, local/source atom and bond mapping, center mapping, and classification
  projection;
- exact classification, wildcard and reversed-bond matching, permissive fallback, precedence,
  missing entries, deterministic diagnostics, and immutable views;
- bundled parameter-set parsing and consistency with registered method requirements.

### C. Built-in method conformance

A manifest-driven suite should cover every built-in method. For each method it should state:

- ID, expected metadata and requirements;
- whether it is parameterless or its test parameter fixture;
- whether it is topology-only or geometry-dependent;
- whether it is reduced-capable and its declared fragment charge policy;
- its basic supported fixture;
- its expected assignment cardinality and charge contract.

The generic suite should validate registration, option-schema validity, applicability, provenance,
finite source-ordered charges, target cardinality, and the relevant topology/geometry behavior.

For all methods other than `formal` and `dummy`, include the neutral-water implementation-sanity
contract described above. If a method cannot meaningfully support water with a minimal test parameter
set, that limitation must be explicit in the manifest and the method must use an equally simple,
chemically interpretable replacement fixture with documented expectations.

### D. Method-specific algorithm tests

Keep a focused file only where it proves logic not covered by generic conformance, for example:

- a method-specific option changes output or iteration behavior;
- a special formal-charge rule is enforced;
- bond parameters, initial charges, or a unique topology are used;
- a fixed small numerical regression is needed to protect a formula.

Numerical values here are implementation regression values. They require explicit fixture,
parameters, options, and tolerance, and should be stable across supported compilers.

### E. Calculation planning and execution contracts

Use synthetic methods where that isolates execution policy from chemistry. Cover:

- deterministic method/parameter ranking;
- full, cutoff, and cover selection, automatic promotion, warnings, and explicit override;
- explicit unsupported mode with no fallback;
- no-plan report state and diagnostics;
- lvalue/rvalue assessment ownership and report lifetime;
- empty collection, empty molecule, singleton target, multiple molecules, and multiple conformers;
- exact source-order cardinality for every execution mode;
- full and reduced validation/solver failures;
- charge correction, fragment target charges, cover ownership, and serial/parallel determinism;
- observer phases, cancellation, and callback isolation.

### F. Adapter and CLI contracts

Adapter tests must prove each format's supported mapping and preservation promise. CLI tests must prove
only application-level behavior: options, exit statuses, file placement, serialized schema, ordered
records, progress rendering, and the declared failure boundary.

## Action plan

The steps below are intentionally incremental. Complete a step with focused tests and a full debug
CTest run before moving to the next. Avoid a broad rewrite of all test files at once.

### Step 0 — Adopt Snitch and make assertions release-safe

**Status: complete.** Snitch is test-only, all C++ tests use active Snitch assertions, and the
release CTest suite executes those checks.

1. Fetch Snitch only from the `CHARGEFW_BUILD_TESTS` branch, following the repository's existing
   find-package-then-`FetchContent` dependency pattern.
2. Convert one focused low-risk target first to named `TEST_CASE`s, including failure output that
   identifies the contract.
3. Verify it in `gcc-debug`, `clang-debug`, `gcc-release`, and `clang-release`.
4. Convert the suite progressively; do not retain standard `assert` checks or local assertion helpers
   whose checks disappear in release builds.

**Exit criterion:** release CTest runs execute active Snitch checks and report the failed named
contract.

### Step 1 — Establish a test inventory and naming convention

**Status: complete.** The A–F inventory is complete, and the oversized calculation, observer,
applicability, and built-in method suites have focused named cases.

1. Classify every existing test as A–F above.
2. Express cases as named Snitch `TEST_CASE`s after observable behavior, not implementation history.
3. Split only the oversized calculation and method executables first; retain existing CTest target names
   where practical to avoid unnecessary build-system churn.
4. Add a short comment at each non-obvious fixture describing the behavior it isolates.

**Exit criterion:** every test has one primary layer and contract; no broad test file contains unrelated
planning, observer, and numerical assertions.

### Step 2 — Complete and remove TODO section 1

**Status: complete.** `test_planning.cpp` owns deterministic ranking, permissive classification,
no-plan diagnostics, and explicit unsupported-policy contracts. Facade/execution tests own the
remaining cardinality, source-order, ownership, and reduced failure contracts; TODO section 1 is closed.

Reconcile the first TODO item with its substantial existing coverage.

1. Move planning selection/ranking cases into a dedicated planning test.
2. Retain one source-order/cardinality test for mixed multiple-molecule/multiple-conformer input per
   execution mode; remove the duplicate result-order assertion from observer coverage.
3. Strengthen explicit unsupported-policy assertions to check the structured execution assessments,
   absent selected plan/policy, and no fallback.
4. Add facade-level singleton input coverage in addition to the existing empty collection and empty
   molecule cases.
5. Keep solver error context with full/reduced execution tests; observer tests should assert only
   notification and cancellation guarantees.
6. Verify lvalue/rvalue assessment ownership using values after the original request has been destroyed.
7. Update `PROJECT.md` if required and remove TODO section 1 only when these cases are focused and
   passing.

**Exit criterion:** section 1 has a small, non-duplicated facade/planning suite and is no longer listed
as unfinished.

### Step 3 — Introduce the built-in method conformance manifest

**Status: complete.** The 22-method manifest owns registry metadata, parameter requirements, option
counts, resource complexity, and reduced-execution capabilities. Its generic workflow covers
applicability, provenance, finite source-ordered charges, topology/geometry fan-out, neutral-water
sanity, atom relabeling, and bond-endpoint invariance. `formal` and `dummy` retain their exact
dedicated contracts; SMP/QEq explicitly uses its minimal water parameter fixture.

1. Convert `test_builtin_methods.cpp` metadata and requirement checks into a table-driven manifest.
2. Add a generic workflow case for all 22 registered methods.
3. Add neutral-water implementation-sanity checks for every method except `formal` and `dummy`.
4. Centralize topology-only versus geometry-dependent fan-out expectations.
5. Centralize atom relabeling and bond-endpoint invariance for deterministic methods.
6. Mark exceptions explicitly in the manifest rather than silently omitting them.

**Exit criterion:** adding a method requires adding exactly one manifest entry plus genuinely
method-specific tests.

### Step 4 — Rationalize method-specific tests

**Status: complete.** Method-specific tests now retain only distinguishing contracts: options,
prerequisites, classification, geometry response, or small explainable numerical regressions. Charge2,
TSEF, VEEM, and MGC have focused water regressions; Charge2 and PEOE verify iteration effects, while
MPEOE verifies bond-attenuation behavior.

1. Preserve only assertions that demonstrate behavior unique to a method.
2. Add a small fixed regression fixture for methods currently lacking one, beginning with `charge2`,
   `tsef`, `veem`, `mgc`, `peoe`, `mpeoe`, `abeem`, and the SQE family.
3. Keep exact values only where the fixture and tolerance are explainable and stable; otherwise assert
   the stronger invariant that captures the intended behavior.
4. Test method options in the owning method test, while generic option-schema validation remains in
   method-options tests.

**Exit criterion:** every built-in method has generic sanity coverage and at least one focused test for
its distinguishing algorithmic behavior.

### Step 5 — Strengthen parameter and feature contracts

1. Expand classifier cases for all supported classification kinds, wildcards, direction-independent
   bonds, precedence, permissive fallback, and failures.
2. Add parameter-view/classification immutability and lifetime tests.
3. Add a bundled-data integrity test over every parameter JSON: parseability, unique IDs, registered
   method ID, and required parameter names.
4. Add prepared-feature cache/lifetime and concurrent read-only access tests as the relevant APIs are
   stabilized.

**Exit criterion:** a parameter matching or bundled-data regression fails at the parameter layer before
it reaches a method calculation.

### Step 6 — Separate reduced-execution mechanics from approximation evidence

1. Create a neutral reduced-execution contract target containing shared cutoff/cover mapping,
   correction, target-charge, ownership, error-context, and serial/parallel tests.
2. Keep whole-molecule-radius equality as an executor regression for all eight supported methods.
3. Label it explicitly as equality under a radius containing the whole molecule.
4. Keep multi-radius error envelopes, disconnected/charged corpus studies, and ChargeFW2 comparison in
   the scientific/compatibility work tracked by `TODO.md`, not as substitutes for executor tests.

**Exit criterion:** reduced executor contracts and approximation-quality evidence are clearly distinct.

### Step 7 — Rationalize adapters and CLI integration

1. Define per-format reader/writer contract tables: source identity, atom order/mapping, conformer
   mapping, malformed-record policy, and preservation behavior.
2. Prefer semantic parse/round-trip checks over output substring checks when a parser is available.
3. Move repeated PDB/mmCIF selection and connectivity expectations into shared contract helpers where
   their public behavior is identical.
4. Add CLI cases as result-schema work lands: invalid request/input, no plan, numerical failure,
   cancellation, output-directory behavior, record ordering, and complete progress-line erasure.
5. Keep parser fuzzing and sanitizer coverage separate as security/robustness gates.

**Exit criterion:** every supported format has explicit, tested import/export promises and CLI tests do
not duplicate library-unit behavior.

## Fixture policy

Use the smallest fixture that demonstrates the contract:

- **unit fixtures:** constructed in code for invalid states and exact topology/geometry;
- **shared molecular fixtures:** water, charged pair, disconnected small systems, and a small asymmetric
  molecule in `tests/cpp/support`;
- **format fixtures:** checked-in files for parser and preservation behavior, with a short description
  of the special syntax or failure they exercise;
- **corpus fixtures:** only where scale or real-format structure is essential.

A fixture must identify whether it is intended for implementation sanity, mapping, numerical
regression, or scientific evidence. Do not reuse water merely because it is available when a charged,
disconnected, or asymmetric case is necessary to expose the behavior.

## Validation cadence for test work

- During edits: build the directly affected `gcc-debug` target and run its focused CTest entry.
- After a coherent change: run `ctest --preset gcc-debug`.
- When changing headers, templates, conversions, or compiler-sensitive tests: run the corresponding
  `clang-debug` target.
- When checks are release-safe or numerical code changes: run focused `gcc-release` and `clang-release`
  tests.
- Use `clang-asan` and `clang-ubsan` for ownership, parsing, mapping, indexing, conversion, and
  numerical-risk changes according to [AGENTS.md](AGENTS.md).

Coverage tools and mutation testing may guide missing-contract discovery, but passing line coverage is
not an acceptance criterion. The measure of success is that each architectural contract has a focused,
active, and non-duplicated test.
