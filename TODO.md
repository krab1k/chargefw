# ChargeFW Delivery Roadmap

This is the actionable delivery tracker for ChargeFW. Technical context and current behavior live
in [PROJECT.md](PROJECT.md); implementation rules live in [AGENTS.md](AGENTS.md); build and test
commands live in [README.md](README.md).

Check an item only after its stated deliverable and validation are complete. Keep atom ordering,
conformer identity, selected method and parameters, options, approximation policy, and diagnostics
explicit at every public calculation boundary.

## Release foundation and public workflow

- [x] Define a high-level autodetect calculation facade that accepts a prepared molecule
  collection, candidate methods, and parameter sets without moving request-specific state into
  `methods::Method`; it selects the highest-priority applicable method and parameter set with
  deterministic ID tie-breaking.
- [ ] Define the corresponding result facade with atom-indexed charges and complete provenance:
  method identifier, parameter-set identity and version, options, conformer, selected policy,
  approximation policy, warnings, and applicability diagnostics.
- [x] Specify deterministic candidate ranking and selection rules, including tie-breaking and the
  behavior when no candidate is applicable; expose an application override rather than silently
  choosing a fallback. The current calculation facade selects the highest-priority applicable
  candidate deterministically and retains applicability diagnostics; caller-directed selection
  remains part of the stable request/result facade work.
- [ ] Add focused API tests for request validation, no-applicable-method outcomes, multiple
  conformers, atom-order preservation, deterministic selection, and provenance completeness.
- [ ] Document the stable request/result boundary and designate it as the sole foundation for
  future CLI, Python, WebAssembly, and MCP adapters.

## Method parity and scientific validation

- [x] Inventory ChargeFW2's `sqe`, `sqeq0`, and `sqeqp` algorithms, required inputs, options,
  failure modes, and all nine associated parameter files without modifying `old/`.
- [x] Implement `sqe` as a stateless method with explicit requirements, prerequisite diagnostics,
  calculation options, and parameter lookup through `ParameterView`.
- [x] Implement `sqeq0` with ChargeFW2-compatible initial-formal-charge behavior and explicit
  diagnostics for unsupported or invalid input.
- [x] Implement `sqeqp` with ChargeFW2-compatible parameterized initial-charge behavior and
  explicit diagnostics for incomplete parameter coverage.
- [x] Migrate the nine archived SQE-family parameter sets to bundled JSON with source attribution,
  schema validation, and stable identifiers.
- [ ] Build a reusable ChargeFW2 comparison fixture that runs old and new methods from identical
  topology, formal charges, conformers, options, and parameter data.
- [ ] Add numerical parity cases for every migrated method and parameter set, including neutral,
  ionic, multi-conformer, and unsupported-parameter cases where applicable.
- [ ] Record per-method numerical tolerances and every intentional behavioral deviation in
  regression fixtures or dedicated compatibility documentation.
- [ ] Verify exact/reference calculations independently from ChargeFW2 `full`, `cutoff`, and
  `cover` modes; do not redesign those modes before compatible behavior and benchmarks exist.
- [ ] Add publication/reference citations and method-specific validation cases for all built-in
  algorithms so scientific assumptions and supported domains are inspectable.

## Parameters, applicability, and diagnostics

- [ ] Define schema-versioning and compatibility rules for parameter JSON if ChargeFW begins
  accepting externally distributed parameter data. Bundled parameter sets remain developer-managed;
  their IDs may be derived from method and name when `metadata.id` is omitted.
- [ ] Add loader tests for malformed JSON, missing required fields, incompatible schema versions,
  duplicate identities, invalid numeric values, and actionable error messages.
- [ ] Test immutable parameter classification across repeated classifications and concurrent
  calculations to prove molecules and cached features are not mutated by matching.
- [ ] Audit every method's `MethodRequirements` and prerequisite checks; add tests that distinguish
  missing geometry, unsupported elements, missing parameters, invalid conformers, and numerical
  limitations.
- [ ] Make parameter-priority selection deterministic and report the full set of considered and
  rejected parameter candidates in calculation diagnostics.

## Core, feature, and result integrity

- [ ] Add invariant tests for molecular graph validation, bond endpoints and orders, conformer
  coordinate cardinality, molecule-collection ordering, and source atom index preservation.
- [ ] Test topology and conformer feature caches for correct invalidation/lifetime behavior,
  multiple conformers, degenerate geometry, and concurrent read-only use.
- [ ] Validate charge results against atom counts, finite values, target total charge, correction
  policy, and conformer identity; retain diagnostic context for failures.
- [ ] Establish numerical robustness tests for singular or ill-conditioned systems, zero-distance
  coordinates, disconnected graphs, and unsupported formal-charge combinations.
- [ ] Benchmark representative small, medium, and large molecular graphs before changing numerical
  solvers, feature representations, or approximation behavior.

## CLI and serializable output

- [ ] Replace the fixed-water demonstration with a user-facing CLI built on the high-level request
  facade, while retaining a focused example or test for the current water workflow.
- [ ] Define CLI inputs for molecule source, method and parameter selection, conformer selection,
  approximation policy, and charge-correction policy; require users to choose or accept a clearly
  reported default.
- [ ] Implement a stable JSON result format containing input record identity, source atom mapping,
  charge assignments, provenance, coverage, warnings, applicability diagnostics, and errors.
- [ ] Add CLI integration tests for successful calculations, no applicable method, malformed input,
  partial parameter coverage, deterministic output, and nonzero exit statuses.
- [ ] Provide machine-readable version information for the library, CLI, parameter data, and JSON
  output schema.

## Integration adapters and batch processing

- [ ] Design the optional RDKit adapter boundary without adding RDKit to `chargefw_core`; cover
  SMILES, SDF/MOL V2000/V3000, Mol2, hydrogen preparation, and conformer transfer policies.
- [ ] Implement the RDKit adapter only after the request/result facade is stable, with explicit
  source atom mapping and no implicit protonation or geometry modifications.
- [ ] Design bounded-memory SDF batch streaming with record-level errors, deterministic output
  order, and a documented policy for continuing after invalid records.
- [ ] Add SDF batch integration tests for ordering, memory bounds, mixed valid/invalid records,
  multi-conformer input, and result-to-source mapping.
- [ ] Design an optional Gemmi adapter with explicit component, model, alternate-location, and bond
  policies for PDB/mmCIF input.
- [ ] Implement and test the Gemmi adapter only when an integration consumer needs it; keep Gemmi
  outside the core library dependency graph.
- [ ] Defer Python bindings until the native request/result and JSON provenance contracts are
  stable; then test ownership, exceptions, atom mapping, and deterministic serialization.

## Packaging, quality, and release engineering

- [ ] Export and install CMake package targets so downstream consumers can use
  `find_package(chargefw CONFIG REQUIRED)` with all required include paths and dependencies.
- [ ] Add a clean downstream consumer configure/build test against the installed package.
- [ ] Ensure installed parameter data is discoverable by installed applications and package
  consumers, with an explicit override for relocatable or custom deployments.
- [ ] Run focused CTest coverage for every change and the complete `gcc-debug` suite before closing
  a scientific-method deliverable.
- [ ] Run the `clang-asan` suite for memory-sensitive, feature-cache, parser/adapter, and numerical
  changes where practical; investigate sanitizer findings before closing the item.
- [ ] Add continuous-integration jobs for GCC debug tests, Clang sanitizer tests, formatting, and
  installed-package smoke tests.
- [ ] Publish versioning, compatibility, and release criteria covering API stability, parameter-data
  compatibility, numerical tolerances, supported methods, and known limitations.

## Product validation and ACC III adoption

- [ ] Establish a reproducible ACC III shadow-comparison corpus with licensed/reference inputs,
  parameter versions, options, expected output mapping, and method-specific tolerances.
- [ ] Compare each migrated method against ACC III and ChargeFW2, classify divergences, and add
  regression tests for accepted behavior.
- [ ] Define acceptance gates for method-by-method ACC III adoption: implementation parity,
  parameter parity, test coverage, provenance, performance, and documented limitations.
- [ ] Conduct independent scientific review of the SQE-family implementation and the numerical
  validation corpus before production adoption.

## Deferred work: start only with a confirmed consumer

- [ ] Redesign cutoff/cover execution only after parity fixtures, performance benchmarks, and an
  explicit reusable execution-policy contract are complete.
- [ ] Add a WebAssembly JSON API only for a browser consumer and only as an adapter over the stable
  native request/result contract.
- [ ] Add an MCP server only as an adapter over stable library or JSON APIs; do not couple it to
  `chargefw_core`.
- [ ] Evaluate ML/GNN estimators only as optional adapters with explicit provenance and coverage;
  do not add them as `chargefw_core` dependencies.
- [ ] Implement custom molecular parsers only where the optional RDKit or Gemmi adapters are
  demonstrably unsuitable.
