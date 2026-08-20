# ChargeFW Delivery Roadmap

This file contains unfinished work only. Implemented state belongs in [PROJECT.md](PROJECT.md), usage
in [README.md](README.md), and implementation rules in [AGENTS.md](AGENTS.md). Remove an item when its
full acceptance criteria are met; do not retain checked history here.

Priorities are ordered by product dependency: stabilize the native contract and reduced execution,
complete the CLI/result boundary, add scientific validation, then package bindings and distributions.

## 1. Calculation and reduced execution

- [ ] Implement deterministic serial cover execution. Define pivot selection, solve halo versus retained
  interior, complete coverage, overlap reconciliation, target charge, correction, diagnostics, and
  method-by-method capability tests before enabling automatic cover.
- [ ] Complete SQE-family reduced validation. Add disconnected and charged SQE+qp fixtures,
  multi-radius full/cutoff convergence/error envelopes, and cover semantics. Treat cutoff/cover as new
  approximations rather than ChargeFW2 parity.
- [ ] Emit explicit-full threshold warnings before fragment/solver allocation as well as retaining them
  in result provenance.
- [ ] Expand facade and planning tests for multiple molecules/conformers, no-plan diagnostics, explicit
  unsupported policies, empty/tiny inputs, fragment solver failures, ownership/lifetimes, and result
  cardinality.
- [ ] Reuse one prepared molecule collection across applicability and calculation in the application
  facade instead of preparing the input twice.
- [ ] Share selected-candidate validation between full and cutoff execution; cutoff currently does not
  reject classifications attached to a method that does not use parameters.

## 2. Result contract, molecular I/O, and CLI

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

- [ ] Build a reusable ChargeFW2 comparison harness using identical topology, formal charges,
  conformers, method options, and parameter data.
- [ ] Reproduce and document ChargeFW2 cutoff behavior for applicable EEM/QEq-like methods, with focused
  fixtures, per-method tolerances, and intentional deviations. Do not force SFKEEM into the generic
  fragment target-charge policy.
- [ ] Add numerical comparison cases for every method and parameter set, including neutral, ionic,
  disconnected, multi-conformer, unsupported-parameter, and failure cases where applicable.
- [ ] Establish a full-versus-cutoff/cover corpus across molecule sizes, methods, radii, charge states,
  and conformers. Report accuracy, charge conservation, runtime, peak memory, and scaling; define
  supported error envelopes rather than implying exact equivalence.
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
