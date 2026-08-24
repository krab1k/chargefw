# ChargeFW Delivery Roadmap

This file contains unfinished work only. Implemented state belongs in [PROJECT.md](PROJECT.md), usage
in [README.md](README.md), and implementation rules in [AGENTS.md](AGENTS.md). Remove an item when its
full acceptance criteria are met; do not retain checked history here.

Priorities are ordered by product dependency: complete implementation contracts, complete the CLI/result
boundary, establish scientific/compatibility evidence, then package bindings and distributions. Numerical
accuracy studies are deliberately separate from implementation-completion work.

## 1. Implementation completion: calculation and reduced execution
- [ ] Expand facade and planning tests for multiple molecules/conformers, no-plan diagnostics, explicit
  unsupported policies, empty/tiny inputs, fragment solver failures, ownership/lifetimes, and result
  cardinality.

## 2. Result contract, molecular I/O, and CLI

- [ ] Record cancellation distinctly at the application result and serialized CLI boundary. A
  cancelled calculation currently has no charges and exit status 1, making its JSON/output behavior
  indistinguishable from a no-plan result. Add an owned cancellation status or diagnostic to the
  stable result schema and provenance, define its exit-status behavior alongside invalid request,
  no-plan, and numerical failure, and test cancellation before a target and during reduced fragment
  work. Do not serialize partial charges.
- [ ] Remove the CLI's full-molecule duplication between `ImportedCollection::molecules` and
  `ImportedExportContext::records[*].molecule`. The current move into assessment depends on retaining
  a second molecule copy for export and makes peak memory roughly twice the imported molecular data
  before prepared features and solver allocations. Redesign imported records/export state to retain
  identity, source mapping, and preservation data without duplicating calculation molecules, while
  preserving source order and generated/preservation-oriented SDF, MOL2, PDB, and mmCIF output. Fold
  this into bounded-memory batch execution and add peak-memory and mapping regression coverage.
- [ ] Make terminal progress rendering erase the complete previous line when switching between target
  and fragment displays or when a shorter message follows a longer one. Keep output on standard error,
  preserve thread-safe rendering, and add a CLI-level assertion that captured progress contains no
  stale suffix characters before the terminal newline.
- [ ] Define the stable result/error schema. Represent successful and failed source records in input
  order with owned import/calculation diagnostics, multi-conformer assignments,
  precision/total-charge rules, and schema compatibility/versioning rules.
- [ ] Make CLI failures structured and distinguish exit statuses for invalid input/request, no
  executable plan, and numerical calculation failure. Do not print a no-plan result only to stdout
  when an output directory was requested.
- [ ] Add bounded-memory batch execution with deterministic record order and memory independent of
  total record count. Preserve the policy that malformed input terminates import unless a separately
  specified recovery mode is introduced.
- [ ] Add complete calculation provenance to SDF, MOL2, and mmCIF in format-appropriate fields. Include
  effective method options once application-facing options exist; do not duplicate JSON schema
  internals unnecessarily.
- [ ] Complete structural mapping/provenance for omitted alternate locations and retain enough opaque
  mmCIF source state for preservation-oriented export.
- [ ] Harden PDB/mmCIF readers with fixtures for empty/incompatible models, unknown elements, insertion
  codes and chain breaks, alternate-location omissions, filtered records, explicit/template/hybrid
  conflicts, and malformed or partially usable multi-block input.
- [ ] Expand parser/CLI integration tests across V2000/V3000 boundaries, malformed/mixed records,
  execution policies, deterministic output, and failure exit statuses.
- [ ] Fuzz native MOL/SDF/MOL2/JSON parsers and run sanitizer coverage before treating arbitrary
  untrusted input as supported.

## 3. Scientific and compatibility validation

These are evidence and release-readiness work, not implementation-completion blockers. Cutoff and cover
remain explicit approximations until their method-specific validation is complete.

- [ ] Establish full-versus-cover validation that proves deterministic pivot and contribution selection,
  complete per-atom ownership, context-only boundary atoms, deterministic overlap reconciliation,
  explicit target-charge/final-correction behavior, and convergence with charge conservation over the
  maintained corpus.
- [ ] Complete SQE-family reduced validation. Add disconnected and charged SQE+qp fixtures,
  multi-radius full/cutoff/cover convergence and error envelopes, and cover semantics. Treat
  cutoff/cover as new approximations rather than ChargeFW2 parity.

- [ ] Build a reusable ChargeFW2 comparison harness using identical topology, formal charges,
  conformers, method options, and parameter data.
- [ ] Reproduce and document ChargeFW2 cutoff behavior for applicable EEM/QEq-like methods, with focused
  fixtures, per-method tolerances, and intentional deviations. Do not force SFKEEM into the generic
  fragment target-charge policy.
- [ ] Validate SMP/QEq archived behavior, fragment target charge, parameters, and numerical convergence
  before enabling any reduced execution; it remains full-only until then.
- [ ] Add numerical comparison cases for every method and parameter set, including neutral, ionic,
  disconnected, multi-conformer, unsupported-parameter, and failure cases where applicable.
- [ ] Establish a full-versus-cutoff/cover corpus across molecule sizes, methods, radii, charge states,
  disconnected systems, and conformers. Report accuracy, charge conservation, runtime, peak memory, and
  scaling; define supported error envelopes rather than implying exact equivalence.
- [ ] Add robustness coverage for singular/ill-conditioned systems, zero-distance and degenerate
  conformers, unsupported formal-charge cases, and fragment solver failures.
- [ ] Audit every method's prerequisites, diagnostics, publication citation, and numerical validation.
- [ ] Add molecular invariant, parameter-loader, immutable-classification, feature-cache lifetime, and
  concurrent read-only tests needed by adapters and parallel execution.

## 4. Python and toolkit integration

- [ ] Define a synchronous toolkit-neutral Python API accepting atomic numbers, formal charges,
  indexed bonds, source names/identities, and zero or more coordinate arrays; return NumPy charge
  arrays plus mappings, provenance, and diagnostics.
- [ ] Add nanobind bindings over the owned facade and test ownership, exception translation, array
  validation, collection ordering, mappings, and all-conformer results.
- [ ] Add a pure-Python `rdkit.Chem.Mol` converter with lazy optional import. Preserve atom indices,
  formal charges, supported bonds, and selected conformers; perform no implicit sanitization,
  hydrogen changes, protonation, embedding, or optimization.
- [ ] Add an explicit helper to attach partial charges/provenance to an RDKit molecule without changing
  formal charges or overwriting properties unless requested.
- [ ] Add Biopython or a native C++ RDKit adapter only after a concrete workflow demonstrates that the
  toolkit-neutral Python route is insufficient. Never make RDKit C++ a core/base-wheel dependency.

## 5. Packaging, installation, and automation

- [ ] Export and install CMake package targets so downstream projects can use
  `find_package(chargefw CONFIG REQUIRED)`; test a clean relocatable consumer build.
- [ ] Add parameter discovery from Python resources, relocatable native installs, and explicit custom
  paths. Normal installed users must not need `CHARGEFW_PARAMETER_DIR`.
- [ ] Build nanobind wheels with bundled parameter data and clean install tests across supported
  CPython/Linux/macOS/Windows targets. The base wheel must not require RDKit.
- [ ] Publish a lazy `rdkit` extra and verify interoperability with supported RDKit distributions.
- [ ] Prepare Conda and versioned OCI image distributions with installation/batch smoke tests.
- [ ] Define machine-readable compatibility versions for library, CLI, parameter data, and JSON schema.
- [ ] Add CI for GCC tests, Clang sanitizers, formatting, native install/consumer tests, wheels, Conda,
  and container smoke tests.

## 6. ACC III adoption gates

- [ ] Establish a reproducible ACC III shadow-comparison corpus with licensed/reference inputs,
  parameter versions, options, mappings, and method-specific tolerances.
- [ ] Classify ChargeFW2/ACC III divergences and retain accepted behavior in regression tests.
- [ ] Define method-by-method release gates for exact parity, approximation accuracy, provenance, I/O
  fidelity, installation, supported platforms, and performance.
- [ ] Obtain independent scientific review of methods and the cutoff/cover validation corpus before
  production adoption.
