# ChargeFW Delivery Roadmap

This file tracks unfinished work in application-delivery order. Current architecture, implemented
methods, and completed milestones belong in [PROJECT.md](PROJECT.md); implementation rules belong in
[AGENTS.md](AGENTS.md); build and test commands belong in [README.md](README.md).

The product priorities are scalable calculation, real-molecule I/O, modern adapters, and easy
distribution. Scientific validation, provenance, diagnostics, and reproducibility are required
acceptance criteria for each priority, not optional follow-up work.

## 1. Calculation contract and scalable execution

- [ ] Complete the application result contract with atom-indexed charges and owned provenance:
  method and parameter-set identity/version, effective options, calculation targets, execution and
  correction policies, radius, warnings, and applicability diagnostics.
- [ ] Add explicit caller selection of method, parameter set, method options, and execution policy.
  Reject unavailable or unsupported selections rather than silently falling back.
- [x] Provide a binding-friendly owned request/result facade that does not expose `Method*`,
  `std::span`, prepared-feature lifetimes, or parameter-storage lifetimes to adapters. Keep the
  existing low-level facade available to advanced native callers; application-facing method options
  remain deferred pending an explicit policy.
- [ ] Introduce a reusable execution policy with explicit `full`, `cutoff(radius)`, and
  `cover(radius)` modes, method capability checks, and no implicit atom-count-based switching.
- [ ] Build shared spatial-neighbor and fragment-mapping support outside `core::Molecule`; preserve
  source atom indices, induced bonds, parameter classifications, and conformer identity in every
  fragment.
- [ ] Reproduce ChargeFW2 cutoff/cover behavior for the EEM/QEq-like methods that previously used
  `EEMethod`, with focused compatibility fixtures and recorded deviations.
- [ ] Define cutoff/cover semantics for SQE, SQE+q0, and SQE+qp, including fragment target charge,
  initial-charge handling, overlap reconciliation, and final total-charge correction. Treat these
  as new approximation policies and validate them against `full`, not as ChargeFW2 parity.
- [ ] Add deterministic parallel scheduling for fragment calculations without global mutable state
  or nested solver oversubscription.
- [ ] Validate result cardinality, finite values, total charge, atom order, molecule/conformer
  identity, mapping, and policy provenance for full and approximate calculations.
- [ ] Add facade and execution-policy tests for invalid requests, explicit selection,
  no-applicable-method outcomes, multiple molecules/conformers, unsupported policy combinations,
  tiny/empty fragments, calculation failures, result ownership, and complete diagnostics.

## 2. Molecular I/O and C++ CLI

- [x] Define the shared import/export contract for record identity, source-to-native atom mapping,
  conformer mapping, diagnostics, unsupported chemistry, and explicit transformations. File-format
  and toolkit types remain outside `chargefw_core`.
- [x] Implement a bounded-memory native MOL/SDF stream adapter for the standalone CLI. Define the
  supported V2000/V3000 subset, reject unsupported/query constructs clearly, preserve source order,
  and report malformed records independently.
- [x] Add a versioned native JSON molecule input adapter for the CLI. Preserve atom, bond, and
  conformer ordering; require atomic numbers and formal charges; do not support atom names or infer
  connectivity from coordinates.
- [ ] Add SDF output that attaches partial charges and calculation provenance without changing
  formal charges, topology, coordinates, or source identifiers.
- [ ] Replace the fixed-water demonstration with a user-facing C++ CLI over the application facade;
  retain the water workflow as a focused example or test.
- [ ] Expose molecule source, method/parameter selection, conformer selection, full/cutoff/cover
  policy, radius, and correction policy through explicit CLI options with reported defaults.
- [ ] Add bounded batch execution with deterministic record order, configurable continue-on-error
  behavior, and memory use independent of total input record count.
- [ ] Implement a stable JSON result format containing input identity, source mapping, assignments,
  full provenance, approximation diagnostics, coverage, warnings, and errors.
- [ ] Add parser and CLI integration tests for valid and malformed records, V2000/V3000 boundaries,
  mixed-success batches, all execution policies, deterministic output, and failure exit statuses.
- [ ] Fuzz native molecular parsers and run sanitizer coverage before treating untrusted input as
  supported.

## 3. Python and toolkit adapters

- [ ] Define a synchronous toolkit-neutral Python API accepting atomic numbers, formal charges,
  indexed bonds, names/source identities, and zero or more coordinate arrays. Return NumPy charge
  arrays plus native mapping, provenance, and diagnostics.
- [ ] Add nanobind bindings over the stable application facade and test ownership, exception
  translation, array validation, collection ordering, atom mapping, and all-conformer results.
- [ ] Add the first toolkit bridge in pure Python for an ordinary `rdkit.Chem.Mol`. Preserve atom
  indices, formal charges, supported bonds, and requested conformers; reject unsupported/query
  chemistry clearly; perform no implicit sanitization, hydrogen changes, protonation, embedding, or
  optimization.
- [ ] Add an explicit helper for attaching calculated partial charges and provenance to an RDKit
  molecule. Validate mapping and conformer identity, never alter formal charges, and do not
  overwrite existing properties unless requested.
- [ ] Keep RDKit chemistry preparation separate from calculation. Any sanitization, hydrogen
  addition/removal, protonation, conformer generation, or optimization helper must be opt-in and
  provenance-recorded.
- [ ] Add a Biopython convenience adapter once a concrete structural-biology workflow defines
  connectivity, component, model, and alternate-location policies; use optional Gemmi support where
  native PDB/mmCIF handling is required.
- [ ] Add a separately built C++ `RDKit::ROMol` adapter only for a demonstrated native consumer;
  never make RDKit's C++ ABI a dependency of `chargefw_core` or Python wheels.

## 4. Distribution and installation

- [ ] Provide parameter-data discovery from Python package resources, relocatable native installs,
  and explicit custom paths; normal users must not need `CHARGEFW_PARAMETER_DIR`.
- [ ] Build nanobind platform wheels with `pyproject.toml`, scikit-build-core, and cibuildwheel (or
  validated equivalents), bundling the native extension and parameter data. `uv add chargefw` must
  work without RDKit.
- [ ] Publish an optional `rdkit` extra so `uv add "chargefw[rdkit]"` installs the seamless RDKit
  workflow. Import RDKit lazily and provide an actionable installation message only when its API is
  requested.
- [ ] Test clean base and RDKit-extra wheel installations with real molecules across the supported
  CPython and Linux/macOS/Windows matrix.
- [ ] Prepare a conda-forge recipe and verify interoperability with Conda-provided RDKit without
  compiling ChargeFW's Python extension against RDKit's C++ ABI.
- [ ] Export and install CMake package targets so native consumers can use
  `find_package(chargefw CONFIG REQUIRED)`, and test a clean downstream build.
- [ ] Publish a versioned OCI/Docker image containing the standalone CLI and bundled parameter data,
  with a batch smoke test and documented input/output mounts.
- [ ] Provide machine-readable versions for the library, Python package, CLI, parameter data, and
  JSON schema, with documented compatibility rules.
- [ ] Add CI jobs for GCC tests, Clang sanitizers, formatting, wheel builds/install tests, native
  install tests, Conda interoperability, and container smoke tests.

## 5. Scientific and performance validation

- [ ] Build a reusable ChargeFW2 comparison fixture that runs old and new methods from identical
  topology, formal charges, conformers, options, and parameter data.
- [ ] Add numerical comparison cases for every method and parameter set, including neutral, ionic,
  multi-conformer, disconnected, unsupported-parameter, and failure cases where applicable.
- [ ] Record per-method numerical tolerances and every intentional behavioral deviation in
  regression fixtures or compatibility documentation.
- [ ] Establish full-versus-cutoff/cover benchmark corpora across molecule sizes, methods, radii,
  charge states, and conformers. Report accuracy, charge conservation, runtime, peak memory, and
  scaling; define supported error envelopes rather than implying exact equivalence.
- [ ] Verify cutoff/cover convergence toward `full` as radius increases and add regression cases for
  fragment boundaries, sparse/disconnected systems, and cover overlap reconciliation.
- [ ] Add numerical robustness coverage for singular/ill-conditioned systems, zero-distance and
  degenerate conformers, unsupported formal-charge cases, and fragment solver failures.
- [ ] Add molecular invariant, parameter-loader, immutable-classification, feature-cache lifetime,
  and concurrent read-only tests needed by I/O, adapters, and parallel execution.
- [ ] Audit every method's requirements and diagnostics for missing geometry, unsupported elements,
  incomplete parameters, invalid conformers, unsupported execution policies, and numerical limits.
- [ ] Add publication citations and method-specific scientific validation cases for all built-in
  algorithms and new approximation policies.

## 6. Product validation and ACC III adoption

- [ ] Establish a reproducible ACC III shadow-comparison corpus with licensed/reference inputs,
  parameter versions, options, source mapping, and method-specific tolerances.
- [ ] Compare each method against ACC III and ChargeFW2, classify divergences, and retain accepted
  behavior in regression tests.
- [ ] Define release and method-by-method ACC III acceptance gates covering exact-method parity,
  approximation accuracy, provenance, I/O fidelity, installation, supported platforms, performance,
  and documented limitations.
- [ ] Conduct independent scientific review of the methods and cutoff/cover validation corpus before
  production adoption.
