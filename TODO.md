# ChargeFW Delivery Roadmap

This file tracks unfinished work in application-delivery order. Current architecture, implemented
methods, and completed milestones belong in [PROJECT.md](PROJECT.md); implementation rules belong in
[AGENTS.md](AGENTS.md); build and test commands belong in [README.md](README.md).

The product priorities are scalable calculation, real-molecule I/O, modern adapters, and easy
distribution. Scientific validation, provenance, diagnostics, and reproducibility are required
acceptance criteria for each priority, not optional follow-up work.

## 1. Calculation contract and scalable execution

- [ ] Complete the application result/export contract. The native JSON writer is the primary complete
  result format and now retains execution/correction policy, radius, threshold, permissive typing,
  and execution warnings alongside record identity/mapping, selected method and optional
  parameter-set IDs, and charge assignments. Add owned method and parameter-set versions, effective
  options, calculation targets, applicability diagnostics, and record-scoped import/calculation
  failures. Propagate applicable provenance to SDF, MOL2, and mmCIF metadata in a separate scoped
  change rather than coupling it to the JSON contract.
- [ ] Add explicit caller selection of method, parameter set, method options, and execution policy.
  Reject unavailable or unsupported selections rather than silently falling back.
- [x] Provide a binding-friendly owned request/result facade that does not expose `Method*`,
  `std::span`, prepared-feature lifetimes, or parameter-storage lifetimes to adapters. Keep the
  explicit native workflow available to advanced callers: applicability returns candidates with
  resolved classifications, selection chooses one candidate, and execution consumes it without
  reclassification. Application-facing method options remain deferred pending an explicit policy.
- [ ] Introduce a reusable execution policy with explicit `full`, `cutoff(radius)`, and
  `cover(radius)` modes, method capability checks, and no implicit atom-count-based switching.
- [x] Build initial shared serial spatial-fragment and mapping support outside `core::Molecule`.
  `features::SpatialFragment` preserves source atom indices, induced bonds, parameter
  classifications, and one source conformer; the shared serial cutoff executor uses it with the
  existing linear neighbor scan.
- [ ] Reproduce ChargeFW2 cutoff/cover behavior for the EEM/QEq-like methods that previously used
  `EEMethod`, with focused compatibility fixtures and recorded deviations.
- [ ] Define cutoff/cover semantics for SQE, SQE+q0, and SQE+qp, including fragment target charge,
  initial-charge handling, overlap reconciliation, and final total-charge correction. Treat these
  as new approximation policies and validate them against `full`, not as ChargeFW2 parity. Cutoff
  now uses induced bonds, method-specific zero/formal initial-charge targets, center retention, and
  explicit final correction with neutral and charged whole-fragment convergence tests; cover overlap
  semantics and multi-radius validation remain open.
- [ ] Add deterministic parallel scheduling for fragment calculations without global mutable state
  or nested solver oversubscription.
- [ ] Validate result cardinality, finite values, total charge, atom order, molecule/conformer
  identity, mapping, and policy provenance for full and approximate calculations.
- [ ] Add facade and execution-policy tests for invalid requests, explicit selection,
  no-applicable-method outcomes, multiple molecules/conformers, unsupported policy combinations,
  tiny/empty fragments, calculation failures, result ownership, and complete diagnostics.

## 2. Molecular I/O and C++ CLI

- [x] Define the shared import/export foundation for record identity, source-to-native atom mapping,
  conformer mapping, diagnostics, unsupported chemistry, and explicit transformations. File-format
  and toolkit types remain outside `chargefw_core`. `ChargeResultDocument` is the initial
  format-neutral result envelope; complete failure/provenance coverage remains tracked in section 1.
- [x] Implement a bounded-memory native MOL/SDF stream adapter for the standalone CLI. Define the
  supported V2000/V3000 subset, reject unsupported/query constructs clearly, preserve source order,
  and report malformed records independently.
- [x] Add a versioned native JSON molecule input adapter for the CLI. Preserve atom, bond, and
  conformer ordering; require atomic numbers and formal charges; do not support atom names or infer
  connectivity from coordinates.
- [x] Add an initial versioned native JSON result writer over the format-neutral result envelope.
  It writes record identity/mapping, selected method and optional parameter-set IDs, assignments,
  and diagnostics; partial-charge values are rounded to at most four decimal places.
- [x] Establish directional native adapter names: `json_input`, `json_output`, `mol_input`,
  `sdf_input`, and `mol2_input`. Future format writers must consume the shared export envelope,
  not serialize calculation results in a front end.
- [x] Add preservation-oriented MOL2 output that copies source bytes and replaces or adds selected
  ATOM partial-charge fields without reconstructing unrelated sections, plus basic generated MOL2
  output from native graph/conformer data when the source format is different. Generated output does
  not infer Tripos atom types or substructure semantics.
- [x] Add preservation-oriented SDF charge output with replace and append modes. Atom-order charge
  vectors are stored in numbered `CHARGEFW_CHARGES_<type-id>` properties without changing formal
  charges, topology, coordinates, source identifiers, or unrelated data fields. Complete calculation
  provenance remains part of the application result/export contract in section 1. Generated output
  supports explicit V2000 or V3000 selection while preservation mode retains the source MOL version.
- [x] Add a Gemmi mmCIF writer. Semantically retain parsed CIF categories and append or replace the
  SB-NCBR partial-charge metadata and charge categories inside selected coordinate-bearing blocks;
  convert PDB through Gemmi and generate one local `UNL` block per nonstructural record. Validate
  atom/model mapping and generated-output round trips. Gemmi serialization is not byte-preserving.
  Future PQR output should be tracked when its adapter is scoped.
- [ ] Replace the molecular-file demonstration with a user-facing C++ CLI over the application
  facade; retain the water workflow as a focused example or test. The current demo accepts
  MOL/SDF/MOL2/PDB/mmCIF/ChargeFW JSON input. All inputs produce JSON and mmCIF; PDB/mmCIF do not
  produce SDF/MOL2
  and expose record-selection and connectivity options. The CLI now has explicit `calculate`,
  `inspect`, `applicability`, `methods`, and `parameters [method-id]` subcommands, including
  method/parameter selection and applicability diagnostics; batch controls, conformer selection,
  complete provenance reporting, and JSON multi-conformer molecular output remain unfinished.
- [ ] Expose molecule source, method/parameter selection, conformer selection, full/cutoff/cover
  policy, radius, and correction policy through explicit CLI options with reported defaults.
- [ ] Add bounded batch execution with deterministic record order and memory use independent of
  total input record count. Reject the entire input on the first malformed record; calculation
  failures remain distinct from input validation failures.
- [ ] Define and complete the stable JSON result schema: ordered entries for successful and failed
  source records; owned import/calculation diagnostics; complete provenance; mapping semantics;
  multi-conformer assignment semantics; precision/total-charge rules; and compatibility/versioning
  rules. Add approximation diagnostics and coverage only with the corresponding execution policies.
- [ ] Add parser and CLI integration tests for valid and malformed records, V2000/V3000 boundaries,
  mixed-success batches, all execution policies, deterministic output, and failure exit statuses.
- [ ] Complete structural-import provenance and mapping. Report selection, altloc, and bond strategy;
  represent omitted alternate locations in source-to-native mapping and diagnostics; and retain
  enough opaque mmCIF source state for preservation-oriented export.
- [ ] Harden Gemmi structural readers with fixture-backed PDB/mmCIF coverage for empty inputs and
  models, incompatible model atom sequences, unknown elements, insertion codes and chain breaks,
  alternate-location omissions, filtered records, explicit/template/hybrid bond conflicts, and
  malformed or partially usable multi-block mmCIF input.
- [x] Keep Gemmi as a required `chargefw_core` dependency. PDB/mmCIF support and its public
  headers, including the Gemmi adapter umbrella, are installed with the core library.
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
- [ ] Add a Biopython convenience adapter once a concrete structural-biology workflow defines how it
  maps onto the existing Gemmi reader's connectivity, component, model, and alternate-location
  policies without duplicating native calculation behavior.
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
  scaling; define supported error envelopes rather than implying exact equivalence. Use
  `tests/fixtures/corpus/cif/10aw.cif` as an initial large structural development fixture; its
  preliminary EEM/QEq full/cutoff reference is recorded in `PROJECT.md` and is not a benchmark.
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
