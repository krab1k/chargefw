# Fixed-Charge Embedding Implementation Plan

Status: EEM implementation authorized, progressing through smaller independently reviewable slices.
Slices 2a-2k are reviewed and complete (2j: `a8f2852`, 2k: `98d8e2e`). Slice 3a SQE+qp/full is
committed as `23b6918`; slice 3b's independent active-charge policy is reviewed and committed as
`5c4b228`. EEM and all three SQE-family methods support full, cutoff, and cover embedding; step 5 is
Astra-reviewed and complete (implementation commit `3000907`). Python slice 6a is Astra-reviewed and
committed as `eb55b97`; newly required steps 6b-6d are not started. CLI step 7 is not started. Embedding support
is limited to EEM and SQE/SQE+q0/SQE+qp by user direction; no other method extensions are planned.
Named-component selection and complete fixed molecular components are required follow-up deliverables
(steps 6b-6d), not implemented by the completed indexed-source slices.
Branch: `fragments`. Research baseline: `4324b11` ([METALS.md](METALS.md)).

## Goal and Scope

Calculate active-molecule charges in the electrostatic field of explicitly prescribed, fixed in-target
ion charges, without requiring fitted parameters for those ions. Removing ions and appending their
charges after an unperturbed calculation does not meet this goal.

The deliverable covers only EEM and SQE/SQE+q0/SQE+qp, with supported cutoff/cover execution and
matching native, Python, and CLI behavior. Start with isolated monatomic sources. Keep parameter-covered
ligands active unless explicitly selected for freezing. Deliver component-name selection for all matching
isolated component instances in imported structures, including monatomic ions such as MG and molecular
components such as SO4 with explicit per-atom charge templates. These are required deliverables, not
optional follow-ons. This scope restriction applies to embedding, not ordinary calculations with other
existing methods.

No new fitted parameter sets, automatic protonation, covalent cutting/capping, fragment-charge generator,
preset library, independent external charge clouds, or response variants for every method.

This file is the user-requested implementation progress tracker. Keep research rationale in METALS.md
and implemented contracts in the owning `docs/` files. Do not duplicate this checklist in TODO.md.

## Working Rules

- **GPT-6 Astra owns architecture, public API design, scientific/coupling decisions, and code review.**
  Astra resolves critical choices before delegating a coding slice; an agent review does not substitute
  for user approval where this plan requires it.
- **GPT-6 Luna owns coding and routine validation** within Astra's agreed design, including documentation
  updates, compilation, tests, formatting, linting, and sanitizer runs. Astra delegates these commands
  to Luna rather than running them itself. Luna reports architectural or scientific ambiguities to Astra
  rather than choosing a new contract or expanding scope.
- For each slice, Astra defines the boundary, Luna implements and runs the relevant checks, and Astra
  reviews correctness, scientific conventions, ownership/mapping, and regressions. Luna addresses review
  findings; Astra verifies fixes by inspecting code and Luna's validation evidence, delegating additional
  checks to Luna. Report the slice ready for user review and stop instead of starting the next step.
- Split numbered steps into smaller review units, normally one commit each. Stop after each slice for
  user review, including API preparation before end-to-end execution. Keep intermediate states usable:
  unsupported combinations must reject explicitly rather than silently omit the field.
- Mark a step complete only after its implementation, focused checks, and review are complete. Record
  its commit and a short validation result below the step; do not accumulate detailed test logs here.
- Follow [AGENTS.md](AGENTS.md) and the commands in [DEVELOPMENT.md](DEVELOPMENT.md). Reuse existing
  solvers, assessment, mapping, and tests; do not build a source-provider framework for hypothetical use.
- Embedding is a scientific approximation independent of `full`, `cutoff`, and `cover`. Default-off
  calculations keep their existing behavior. Automatic selection must compare the same partition.

## 1. Settle the Small Contract

- [x] User authorized starting with EEM and explicitly requested smaller review boundaries such as API
  preparation. The EEM request API slice below is authorized.

### Proposed C++ Shape

Keep one request-level option, independent of method options and execution policy:

```cpp
struct FixedAtomCharge {
    std::size_t molecule_index; // collection index, zero-based
    std::size_t atom_index;     // atom index within the original molecule, zero-based
    double charge;
};

struct FixedChargeEmbedding {
    std::vector<FixedAtomCharge> sources;
    std::string charge_provenance;
};

struct AssessmentRequest {
    // existing fields...
    std::optional<FixedChargeEmbedding> fixed_charge_embedding;
};
```

For each affected target, calculate `Q_active` as the formal-charge sum of unselected atoms in the
prepared active molecule. Prescribed fixed charges define the source sites but do not alter `Q_active`;
the modeled total is `Q_model = Q_active + sum(fixed charges)`. Preserve the original supplied formal-charge
sum as audit provenance even when it differs from `Q_model`. Do not infer charge from element/component ID,
apply CCD rules, modify input atoms, or change prepared active formal charges.

### Contract Decisions Proposed For Approval

| Decision | Proposed contract |
| --- | --- |
| Scope and selection | Sources belong to the same input molecule/target as active atoms; no cross-record or external point cloud. Each source is an explicitly indexed atom with a finite real prescribed charge. Initial support is isolated atoms only: reject duplicate indices, out-of-range indices, selected atoms with any graph bond, and a target with no active atoms. An empty source list normalizes to embedding disabled before capability checks and produces no embedding provenance. No element-based inference. |
| Coordinates and identity | Use original target conformer coordinates for each calculation target. A source index identifies the same atom across conformers; evaluate its position in the corresponding conformer. Preserve original molecule, conformer, and atom order in returned results. |
| Charge model | `Q_active` is the formal-charge sum of unselected atoms in the prepared active molecule; selected atoms' imported formal charges do not affect it. The modeled total is `Q_active + sum(fixed charges)`, and fixed values appear exactly once in reassembled output. Preserve the original supplied formal-charge sum as audit provenance, not as a constraint on the modeled total. Do not alter prepared input atoms. |
| Validation | Reject nonfinite charges/totals/coordinates, coincident active/source sites, and invalid graph/source scope. Do not clamp distances. Do not define a universal near-contact radius or warning absent a method-specific scientific basis. |
| Ownership | `AssessmentRequest` owns raw selectors. `AssessmentResult` owns its validated partition, source values, mappings, and geometry lifetime alongside its current molecule/prepared-feature owners. Its facade execution passes internally validated target-local context tied to the active prepared data and candidate. `CalculationInput` receives a read-only non-owning `std::span<const FixedPointSource>` (each source has `core::Position position` and `double charge`), valid throughout the method call. Coupling stays method-specific. |
| Lower-level execution | Keep the existing public `CalculationRequest` unchanged and non-embedding initially: it consumes already prepared/classified data and cannot accept raw selectors or repartition. Only the owned assessment facade partitions, validates, executes with its tied internal context, and reassembles results at the facade boundary. |
| Capability/resource | Add an explicit embedding capability to method requirements/assessment, independent of `full`/`cutoff`/`cover`; only EEM and the SQE family implemented in steps 2-3 advertise it, and only full execution until step 5. Unsupported candidates/modes are rejected, never silently downgraded. Resource thresholds continue to apply per active solve molecule under the existing policy; do not add a source-work estimate or count override. |
| Provenance | Add structured embedding provenance to effective/result output: source molecule/atom indices and exact prescribed charges, caller-supplied charge-provenance label, and per-target original and active totals. The selected method, parameter set, and options already identify the coupling; do not add a separate kernel field. JSON may summarize imported structural component labels without assigning components. Keep provenance separate from execution mode/radius. The existing native `ExecutionResult` and adapter `ChargeCalculationResult` are the result boundaries; JSON serializes the structure, molecular charge arrays remain original-order. |

`CalculationInput` is a public method-level API, so adding its source span is a public API change, not
an internal-only detail. It carries numerical inputs, not original-atom selection or result reassembly.
Callers invoking `Method::calculate()` directly must satisfy the method's declared requirements,
including embedding capability; this is a low-level API precondition, not a guarded dispatch interface.
The facade rejects unsupported embedding candidates before classification/planning and checks the
selected method/mode in its existing execution validation when source plumbing is connected. It must
never discard the environment to obtain a usable candidate. Supporting methods validate numerical
source data they consume; unrelated algorithms do not receive embedding guards. No new dispatch wrapper
is needed. The facade also owns partition, budget, and mapping validation. Source buffers must be owned
by the assessment or the individual execution, never shared mutable scratch state, so concurrent plan
reuse remains safe.

No CCD charge assignment, fragment-charge generation, automatic chemistry policy, or core molecule
change is part of this proposal. Python should expose the same zero-based `(molecule_index, atom_index,
charge)` values and charge provenance through `assess`/`calculate` and plan snapshots; CLI `calculate`
and `applicability` should accept a simple repeatable source selector with the same explicit indices,
charge, and provenance. Both use prepared active formal-charge sums; original input sums remain audit
provenance. Exact Python/CLI spelling remains
for their implementation steps and is not a step-1 selector-syntax approval; neither may parse chemistry
or create its own partition policy.

Review outcome: user authorized starting EEM implementation, with API preparation as a smaller review
boundary. Implement the contract incrementally without enabling incomplete execution paths.
Astra adopted the independent active-charge model for the
current 3b follow-up; see the dated evidence note in METALS.md.

The initial contract above describes the delivered indexed, isolated-atom interface. Steps 6b-6d extend
source preparation to named component instances and whole fixed molecular components. Their public API
and selector syntax require review before implementation; no automatic charge assignment is authorized.

## 2. Deliver EEM Full Execution End to End [x]

- [x] Implement one complete native path from assessment through reconstructed output in slice 2k.

### Review Slices

- [x] **2a: Native request API preparation.** Add `FixedAtomCharge`, `FixedChargeEmbedding`, and the optional
  `AssessmentRequest::fixed_charge_embedding` field. Empty source lists normalize to disabled; reject every
  nonempty embedding request explicitly before preparation/classification until execution is supported.
  Test both assessment ownership overloads, default/empty compatibility, and explicit/automatic
  selection rejection. Document this safe API boundary. Do not add numerical method input, partitioning,
  capabilities, provenance, bindings, or solver changes in this slice.
- Slice 2a validation: GCC debug focused planning test and full suite (56/56), Clang debug focused
  planning test, affected clang-tidy target, and whitespace checks passed. Luna implemented; Astra
  reviewed and verified the strengthened rejection tests and full GCC suite. Commit: `8749313`.
- [x] **2b: Source-selection validation.** Validate original molecule/atom index bounds, duplicate
  selections, finite prescribed charges, isolated selected atoms, and at least one remaining active
  atom per affected molecule. Report invalid selections before the existing unsupported-execution gate;
  valid nonempty selections still reject there. Accept finite zero/fractional charges without element
  inference or formal-charge mutation. Keep empty/default requests unchanged. Geometry, budget checks,
  owned partitioning, method capability, and solver execution remain subsequent work.
- Slice 2b validation: GCC debug full suite (56/56); focused GCC/Clang debug, ASan, and UBSan planning
  tests; affected clang-tidy target; formatting and whitespace checks passed. Luna implemented; Astra
  reviewed and verified final GCC regressions, including both bond endpoints. Commit: `9d787b9`.
- [x] **2c: Geometry and charge-budget validation (original total policy).** For each affected molecule, require conformers and
  finite coordinates for every atom in every conformer; reject exact active/source coincidence without
  a distance clamp or near-contact threshold. Validate finite prescribed-source sums and the budget
  `Q_active = Q_original - sum(fixed charges)`, using the original supplied formal-charge total. Keep
  source selection, geometry, and totals target-local; do not infer formal charges or require prescribed
  charges to equal them. This original-minus-source budget rule is superseded by the independent active
  charge policy adopted in the 2026-10-04 follow-up; the validation record below describes the 2c
  implementation at that time. Valid nonempty requests still reach the unsupported-execution gate. No owned
  partition, solver, method capability, or public API change belongs to this slice.
- Slice 2c validation: GCC debug full suite (56/56); focused Clang debug, GCC/Clang release, ASan, and
  UBSan planning tests; affected clang-tidy target; formatting and whitespace checks passed. Luna
  implemented; Astra reviewed and verified final GCC regressions, including conformer isolation.
  Commit: `6050b20`.
- [x] **2d: Private owned partition construction.** Add a private factory that owns an active collection
  with original molecule/conformer ordering, active-to-original atom/bond maps, exact fixed-source
  charges and per-conformer positions, original/active budgets, and the charge-provenance label.
  Preserve atom metadata and formal charges; retain disconnected active components together. Copy
  unaffected molecules with identity mappings. Share existing validation with assessment, without
  constructing a discarded partition on its still-unsupported path. Test ownership after input
  destruction/moves and ordering/bond/geometry/budget mappings through the private interface. Do not
  attach the partition to assessment owners or enable classification, reconstruction, or solvers yet.
- Slice 2d validation: GCC debug full suite (57/57); focused Clang debug and GCC/Clang release tests;
  ASan/UBSan calculation suites (9/9 each) and final focused reruns (2/2 each); affected clang-tidy
  targets and whitespace checks passed. Luna implemented; Astra reviewed ownership and mappings and
  verified final regressions. Commit: `219628d`.
- [x] **2e: Method input and capability declarations.** Add `methods::FixedPointSource` (position and
  prescribed charge) and a trailing default-empty borrowed source span to `CalculationInput`, preserving
  its existing temporary-owner protections. Add a default-false scientific capability on
  `MethodRequirements`, separate from execution resources. Keep all built-ins unsupported, including EEM.
  Leave all algorithm implementations unchanged. Direct callers must check declared capability before
  supplying nonempty sources. Test borrowed source lifetime/identity semantics and unsupported capability
  declarations across the registry, not rejection of direct calls that violate preconditions. The
  current facade guard continues rejecting nonempty requests. Candidate/mode checks belong in existing
  assessment and execution validation when embedding is connected, not in new wrappers. No solver,
  assessment ownership, facade source plumbing, or reduced-execution support is enabled here.
- Slice 2e revision: removed per-algorithm guards and their helper at user request; retained source
  input, capability flag, and lifetime/declaration tests. All built-in implementations remain unchanged.
  GCC debug full suite (57/57), focused GCC/Clang debug tests (3/3 each), sequential ASan/UBSan focused
  tests (2/2 each), affected clang-tidy targets, and whitespace checks passed. Luna implemented; Astra
  reviewed the reduced scope and verified final regressions. Commit: `5f2c935`.
- [x] **2f: EEM numerical fixed-source response.** Enable EEM's numerical capability and subtract
  `sum(kappa * source_charge / distance)` from its atomic RHS, preserving the active-active matrix and
  caller-supplied active target budget. Use the selected set's `kappa`; require no source parameters.
  Validate basic numerical source inputs at the boundary. Use existing `core::distance` and ordinary
  EEM arithmetic without extra per-pair or intermediate-overflow checks. Direct callers must supply
  noncoincident source/active positions; the facade retains its existing coincidence validation.
  Keep empty-source behavior unchanged.
  Test an independent two-atom reference, nonunit scaling, source sign/position and budget behavior,
  zero sources, invalid inputs, and unchanged source storage. Only EEM consumes sources; leave unrelated
  algorithms unchanged. Assessment still rejects nonempty embedding requests until the remaining facade
  ownership, capability/mode checks, provenance, and reconstruction are connected.
- Slice 2f structure after user review: basic numerical source/geometry checks live in a private shared
  validator called only by supporting methods. A small EEM-local helper sums the point potential;
  `calculate()` remains focused on assembly, RHS adjustment, and solving. Removed special `hypot`,
  per-operation checks, and extreme-arithmetic tests to stay consistent with existing EEM. No capability
  guard, dispatch wrapper, kernel framework, or duplicate pair-distance pass was added.
- Slice 2f revised validation: GCC debug full suite (57/57); focused GCC/Clang debug and release tests
  (3/3 each); sequential ASan/UBSan focused tests (3/3 each); affected clang-tidy targets and whitespace
  checks passed. Luna implemented; Astra reviewed the simplified kernel and boundary checks.
  Commit: `0e28b30`.
- [x] **2g: Active-partition assessment.** Own the private partition at a stable address in assessment,
  prepare its active collection, and retain original molecules for public identity. Reject methods
  lacking embedding capability before parameter classification. Assess capable methods on the active
  graph with unchanged matching, map structured diagnostic indices back to original atoms/bonds, and
  make any active-index explanatory text explicit. Return structured unsupported-execution rejections
  for otherwise-applicable embedded candidates in every mode; create no embedding execution plans yet.
  Apply the temporary block in the existing per-mode planning check, without overwriting candidate
  execution assessments or enumerating execution modes a second time.
  Test missing-ion-parameter exclusion, strict active coverage, original mappings, ownership/moves,
  explicit and automatic selection, and empty/default compatibility. No solver dispatch, provenance,
  or result reconstruction is connected in this slice.
- Slice 2g validation: full GCC debug suite (57/57), focused Clang debug and GCC/Clang release tests
  (3/3 each), sequential ASan/UBSan calculation suites (11/11 each), affected clang-tidy targets, Ruff,
  and whitespace checks passed. Luna implemented; Astra reviewed ownership, capability filtering, and
  diagnostic mapping. Final integration cleanup passed focused tests and full GCC. Commit: `491003e`.
- [x] **2h: Private result reconstruction.** Add a small private helper using the validated partition
  to scatter active charges to original atom indices and insert exact fixed values. Preserve assignment
  order, molecule/conformer targets, and method/parameter identity. Trust factory-owned mapping
  invariants; check active charge-vector dimensions and use bounds-checked indexing rather than adding
  another partition validation pass. Test multiple targets/conformers, interleaved sources, identity
  partitions, and unchanged totals/source values. Do not connect execution or public provenance yet.
- Slice 2h validation: full GCC debug suite (57/57), focused GCC and Clang debug/release tests,
  sequential ASan/UBSan calculation suites (9/9 each), affected clang-tidy targets, and whitespace
  checks passed. Luna implemented; Astra reviewed scatter mappings, exact fixed-value insertion, and
  metadata preservation and verified full GCC regressions. Commit: `255e0e2`.
- [x] **2i: Native and JSON embedding provenance.** Add optional structured embedding metadata to
  `EffectiveCalculation`: original-index sources and charges, caller provenance label, interaction-model
  identifier, and original/active charge totals once per original molecule. Serialize it under effective provenance
  only when present; ordinary JSON remains unchanged. Extend the result schema with the optional field
  and test synthetic successful/failed results without enabling plans. Adapter checks cover structural
  metadata validity, not a second scientific assessment. Keep Python typed snapshots/request options,
  CLI source input, and execution integration for later slices.
- Slice 2i validation: full GCC debug suite (57/57), focused GCC/Clang debug JSON tests, sequential
  ASan/UBSan JSON tests, affected clang-tidy target, and whitespace checks passed. Emitted synthetic JSON
  validated against the updated schema using installed jsonschema. Luna implemented; Astra reviewed
  optional-field behavior and corrected test budget/source consistency. User requested renaming `kernel`
  to `interaction_model` and `budget` to `charge_totals` with explicit original/active total field names.
  Focused GCC JSON output test (1/1), full GCC suite (57/57), schema meta-validation, and whitespace
  checks passed after terminology update. Commit: `2e5df43`.
- [x] **2j: Private full-execution source plumbing.** Extend only the private full-execution entry point
  with a trailing optional partition pointer. Require the factory-owned active molecule collection to
  match the prepared molecules by object identity, and check the selected method's embedding capability
  at full-execution validation. For each target, pass its active charge budget and materialize borrowed
  point sources from the matching conformer's validated partition positions for the synchronous method
  call. Keep ordinary no-partition calls unchanged; leave facade planning, result reconstruction,
  provenance, and reduced execution disconnected. Test budget/source target mapping, conformer identity,
  unsupported methods and mismatched ownership. Astra reviewed and accepted the implementation.
- Slice 2j implementation validation: GCC debug full suite (57/57); focused `test_fixed_charge_partition`
  passed under Clang debug and GCC/Clang release (1/1 each); sequential ASan/UBSan runs of
  `test_calculation`, `test_planning`, `test_fixed_charge_partition`, `test_observer`,
  `test_execution_policy`, `test_calculation_targets`, `test_reduced_execution`,
  `test_representative_execution`, and `test_cover_execution` passed (9/9 each). Affected clang-tidy
  targets, formatting, and whitespace checks passed. Astra reviewed the fixes, including serial-only
  capturing-method checks, the independent parameterized EEM reference, and source lifetime/identity
  handling. At the 2j boundary, facade tests verified nonempty embedding requests produced no plans.
  Commit: `a8f2852`.
- The resource-accounting detour was dropped after user review; 2k integrates native EEM full-mode
  facade execution. 2k review removed the redundant interaction-model result field introduced in 2i;
  method ID/options remain the coupling provenance.
- [x] **2k: Native EEM full-mode facade and JSON component summaries.** Allow only full plans for
  validated fixed-charge EEM requests. Reuse one private dispatch/emitter path, produce embedding
  provenance before execution, and reassemble successful results once at the facade boundary. Preserve
  unsupported reduced modes, automatic warning policy, ordinary direct requests, and metadata on
  failure/cancellation. JSON retains flat sources as authority and may group consistently labeled
  monatomic sources by component ID and exact prescribed charge; labels come from imported structural
  metadata, never CCD assignment. Do not expand Python/CLI request APIs. Astra reviewed and accepted the
  implementation; commit recorded in git history.
- Slice 2k validation: GCC debug full suite (57/57). Native facade suites (`test_builtin_methods`,
  `test_planning`, `test_calculation`, `test_observer`, `test_json_output`) passed GCC/Clang debug and
  GCC/Clang release (5/5 each); sequential ASan/UBSan with `test_fixed_charge_partition` passed (6/6
  each). Component adapter suites (`test_json_output`, `test_mmcif`, `test_mmcif_output`, `test_pdb`)
  passed Clang debug and GCC/Clang release and sequential ASan/UBSan (4/4 each). Component schema
  instances validated. Affected clang-tidy targets completed with an existing Gemmi `MmcifReader`
  special-member warning; formatting and whitespace checks passed.
- Named-component selection is now required in steps 6b-6d. CCD-backed charge assignment remains separate
  and is not authorized by component selection alone.
- Step 3, the SQE family, remains the next implementation scope after user review.

Primary files: `include/chargefw/calculation/{assessment,calculation}.h`,
`src/calculation/{assessment,calculation,full_execution,target_execution}.*`,
`include/chargefw/methods/{calculation_input,method_requirements}.h`, and
`src/methods/builtin/eem.cpp`, `src/adapters/native/json_output.cpp`,
`schemas/result-1.0.schema.json`, and the calculation/JSON adapter tests.

- Prepare the partition once before candidate checks. Own active molecules and their features with safe
  lifetimes; retain disconnected active components together and preserve original conformer identity.
- Declare embedding capability explicitly. Initially enable EEM/full only; automatic planning must not
  choose an unsupported method or reduced mode, or override resource policy to manufacture a plan.
- Enforce embedding compatibility in candidate assessment before classification and in existing
  execution validation before dispatch. Direct method callers are responsible for declared preconditions;
  do not add guards to algorithms that do not implement embedding or introduce a separate embedding-only
  dispatcher.
- Assess prerequisites and classify active atoms/bonds only with the unchanged matcher. Translate active
  diagnostic indices back to original atom/bond indices.
- Pass the prepared active formal-charge sum to the solve and subtract
  `sum(kappa * fixed_charge / distance)` from its RHS,
  using the selected EEM set's `kappa`. Do not read ion atom parameters.
- Reassemble original-order charges once at the common result boundary. Fixed values are copied exactly;
  reusable plans, progress, cancellation, and result identity retain their existing guarantees.
- Add structured embedding provenance while retaining the existing resource policy for per-active-molecule
  thresholds; do not add a source-work estimate or resource count override.
- Document the supported native/full slice in `docs/PROJECT.md` and `docs/NATIVE.md`.

Focused check: one small active pair plus a missing-parameter ion, compared with an independently
assembled constrained EEM system. Demonstrate nonzero active response, exact fixed charge, the modeled
total as active plus prescribed charges, and reusable original-order output. Add a compact invalid-request
table to existing planning tests for crossing bonds, invalid selections/geometry, unsupported modes, and
missing active coverage.

## 3. Add the SQE Family

- [x] Complete SQE-family full embedding with its shared source field and independent active-charge policy;
  Astra-reviewed and complete.

Primary files: `src/methods/builtin/{sqe,sqeq0,sqeqp}.cpp` and method prerequisite checks.

- [x] **3a: Shared source field and SQE+qp full capability.** Validate sources before SQE-core empty/no-bond
  returns. Subtract the active-width Gaussian-to-point potential before transfer projection, preserving
  `-H*s + D*s`, active-active terms, and bond hardness. Enable full embedding capability only for SQE+qp;
  keep SQE and SQE+q0 disabled pending budget compatibility. Preserve SQE+qp's once-global seed
  normalization to `Q_active` and component seed totals. Validate scalar RHS, width sign/zero, invalid
  source inputs before no-bond returns, disconnected components, and facade source mapping. Astra reviewed
  field sign and seed preservation, including the strengthened facade scalar reference.
- Slice 3a validation: GCC debug full suite (57/57); focused `test_sqeqp`, `test_sqe`,
  `test_builtin_methods`, `test_planning`, `test_calculation`, and `test_fixed_charge_partition` passed
  under Clang debug and GCC/Clang release (6/6 each) and sequential ASan/UBSan (6/6 each). Affected
  clang-tidy targets, formatting, and whitespace checks passed. Commit: `23b6918`.
- [x] **3b: Independent active-charge policy and SQE/SQE+q0 full capability.** Calculate the active target
  from prepared formal charges on unselected atoms; keep selected imported formal charges in original-total
  audit provenance only, and report modeled total as active plus prescribed sources. Remove obsolete
  prerequisite/runtime budget gates and public budget-selection API; retain SQE neutral-component
  prerequisites, SQE+q0 formal seeds, and SQE+qp global seed normalization. Preserve finite prescribed
  source-value and source/active-geometry checks. Astra reviewed and accepted; complete.
- Slice 3b validation: GCC debug full suite (57/57); focused `test_fixed_charge_partition`,
  `test_calculation`, `test_planning`, `test_sqe`, `test_sqeqp`, `test_builtin_methods`,
  `test_method_prerequisites`, `test_method_applicability`, and `test_json_output` passed under Clang
  debug and GCC/Clang release (9/9 each). The same set plus `test_observer` passed sequential ASan/UBSan
  runs (10/10 each). Affected clang-tidy targets completed with an existing Gemmi `MmcifReader`
  special-member warning; schema meta-validation, formatting, and whitespace checks passed.

Step 4 was removed when embedding scope was restricted to EEM and the SQE family. Later step and slice
numbers remain unchanged to preserve references to completed work.

## 5. Propagate Embedding Through Cutoff and Cover

- [x] Complete step 5 facade wiring and method qualification; Astra-reviewed and complete.
- [x] **5a: Private fragment source-input plumbing (`71495e8`).** Implemented as a default-empty borrowed source span
  on `detail::calculate_fragment_charges`; the complete source list is forwarded unchanged to the method
  call. Focused tests cover EEM and SQE+qp, including a source outside the fragment radius and full-atom
  fragment agreement with direct full calculation. Astra reviewed source forwarding, lifetime, and
  reference-target semantics; accepted and complete. This slice alone did not enable reduced embedding
  through assessment or facade execution; step 5b wires those paths.
- Slice 5a validation: focused `test_reduced_execution` passed under GCC debug, Clang debug, GCC/Clang
  release, ASan, and UBSan; full GCC debug suite passed (57/57). The affected clang-tidy target, format,
  and whitespace checks passed.
- [x] **5b: Executor wiring and method qualification.** Implemented cutoff/cover partition forwarding,
  complete conformer-specific source materialization, active-to-original source-index diagnostics, and
  facade mode planning for EEM and the SQE family. Tests cover remote sources, full/whole-active-radius
  agreement, unaffected targets, multi-conformer mapping, automatic cutoff/cover selection, conservation,
  and cancellation/reuse. Astra reviewed and accepted; complete.
- Slice 5b validation: full GCC debug suite (57/57); focused `test_planning`, `test_observer`,
  `test_reduced_execution`, `test_fixed_charge_partition`, `test_cover_execution`, `test_calculation`,
  `test_representative_execution`, `test_builtin_methods`, `test_sqe`, `test_sqeqp`, and
  `test_electronegativity_equalization` passed under Clang debug and GCC/Clang release, plus sequential
  ASan/UBSan. The final parameterized `test_observer` cancellation suite also passed under GCC debug,
  Clang debug, GCC/Clang release, ASan, and UBSan. Affected clang-tidy targets and format/whitespace
  checks passed.

Primary files: `src/calculation/{reduced_execution,cutoff_execution,cover_execution}.cpp`.

5b qualification requirements:

- Build references and conserved groups on the active graph using the prepared active formal-charge total.
  Preserve the one-time SQE+qp normalization from step 3.
- Supply the complete fixed-source environment to every fragment, including sources beyond its radius.
  Start by reusing direct source evaluation; precompute/project potentials only if necessary, without
  introducing a cache framework or changing coupling conventions.
- Apply conservation corrections to active results before source reinsertion. Never include fixed
  sources as fragment variables, correct their values, or double-count their contribution.
- Preserve source-index diagnostics and finite-radius cover's existing order-dependent pivot behavior.
- Update execution capability/resource handling and `docs/PROJECT.md` for supported combinations.

Focused 5b evidence is in `test_planning` and `test_observer`: the facade compares cutoff/cover with full
execution for EEM and the SQE family, including remote conformer-specific sources and disconnected
SQE+qp component references; observer tests cover fragment cancellation and repeated execution.

## 6. Expose the Policy in Python

- [x] **6a: Python request, binding, and result provenance.** Implemented for assessment, direct calculation,
  and reusable plans; Astra-reviewed and committed as `eb55b97`. This completes the indexed
  Python interface; required named-component work in 6b-6d remains open.

Primary files: `python/chargefw/{calculation,_calculation_options}.py`, `python/src/calculation.cpp`,
associated value types, exports, and extension stubs.

Translate the approved molecule/atom-index request into native data. Keep immutable request snapshots,
plan/collection checks, conformer handling, and adapter source mappings. Python must not partition the
graph or implement its own charge/capability policy. Update `docs/PYTHON.md`.

Focused check: one public workflow exercising assessment, direct calculation, and plan reuse over a
multi-conformer input, with matching native results and provenance. Run existing calculation/output
and typing checks rather than duplicating numerical tests in Python.

Validation: all nine registered Python CTest suites passed with
`ctest --test-dir build/gcc-debug --output-on-failure -E '^cpptest$' -R '^test_chargefw_python_'`
(including `test_chargefw_python_mypy`). The focused calculation test also passed under Clang debug,
GCC/Clang release, ASan, and UBSan; the binding target passed clang-tidy. Ruff checks and format checks
passed for the modified Python files.

### 6b. Select Monatomic Sources by Component Name

- [ ] Add a structural preparation interface selecting every matching component instance by imported
  component ID, for example MG, with an explicit prescribed charge per monatomic instance (`MG -> +2`).
  Resolve names into original molecule/atom indices, then use the existing native embedding path.
- [ ] Define component-instance identity and target/conformer scope using existing import mappings.
  A component ID selects instances, not one combined group of all atoms sharing that ID. Reject matching
  instances with graph bonds to outside atoms; do not silently skip them or delete coordination bonds.
- [ ] Validate monatomic identity, element consistency, finite charges, duplicate/overlapping selections,
  and actionable unmatched-selector diagnostics. Settle exact API and selector syntax during review.
  Imported component names select sources; neither element names nor missing imported formal charges
  silently assign them. CCD-backed assignment would require a separate explicit policy.

### 6c. Support Whole Fixed Molecular Components

- [ ] Extend native source validation to permit internal bonds within a completely selected fixed
  component while rejecting every active/fixed crossing bond. Retain at least one active atom per target.
  Internal fixed bonds do not enter the active solver or require active bond parameters.
- [ ] Require a complete finite per-atom charge template keyed by imported component atom names, with
  declared total, tolerance, and charge-model provenance. Match every atom unambiguously; reject missing,
  extra, or ambiguous template entries and incompatible component identity/chemical state. Do not infer
  a molecular partial-charge distribution from formal charges or the component total.
- [ ] Apply the template to every selected instance, including SO4, and preserve original mappings,
  conformers, exact fixed values, and independent active-charge accounting. Reuse existing EEM/SQE source
  coupling in full, cutoff, and cover; add no new solver or automatic charge generator.

### 6d. Integrate and Validate Named-Component Workflows

- [ ] Expose structural preparation consistently to native and Python callers, feeding the shared native
  assessment request rather than duplicating partition or scientific policy in bindings. Preserve the
  existing indexed-source API for callers with explicit mappings.
- [ ] Test an mmCIF workflow with repeated MG and SO4 instances, explicit ion charges and a caller-supplied
  sulfate template, assessment/calculation/plan reuse, conformers, and original-order output. Verify
  equality with equivalent indexed sources and nonzero active response. Cover crossing bonds, incomplete
  templates, unmatched names, target isolation, and unsupported methods. A test template is not a
  validated sulfate preset.
- [ ] Document the implemented workflow and add an executable example. Include source-assignment
  provenance and resolved original atom mappings; component summaries must not replace exact sources.

Steps 6b-6d are required before final qualification. Deliver named ions first, then molecular templates;
neither requires a preset library, CCD charge inference, external charge clouds, or covalent cutting.

## 7. Expose the Policy in the CLI and Serialized Results

- [ ] Add indexed and component-name source input to both `calculate` and `applicability`, including
  explicit monatomic charges and complete molecular templates, with matching output provenance.

Primary files: `apps/chargefw/{input,output}.cpp`, shared result serializers and their schemas.

Use unambiguous molecule/atom scope and index base for the selector syntax settled during CLI
implementation. Reuse the structural component resolver from steps 6b-6d; keep component-ID selection
simple rather than adding a general selection language. Translate resolved sources into the same native
assessment request and do not infer Mg2+ from PDB/mmCIF element names. Retain normal charge arrays and original
atom mappings in every output path. Update `docs/CLI.md`, `docs/FORMATS.md`, and one concise executable
usage example.

Focused check: installed-CLI round trips for indexed sources and named MG/SO4 components through both
applicability and calculation, plus malformed-input coverage in existing argument tests. Check structured
provenance against the output schema. Use CTest
installation fixtures rather than parameter-dependent commands from the build tree.

## 8. Review and Qualify the Initial Deliverable

- [ ] Complete integration review and record validation evidence and scientific limits.

Reuse the small fixtures above. Establish a rebuilt baseline before implementation; use focused
`gcc-debug` checks during each step and the full debug suite after coherent changes. Follow repository
cadence for Clang header checks, release numerical checks, sequential ASan/UBSan mapping/lifetime
checks, and clang-tidy. Run the full required matrix at the final integration milestone, not after
every small edit. Format changed C++ files with clang-format throughout.

Run one reproducible small ion/ligand distance scan with explicitly prepared chemistry, comparing
embedded and ion-omitted results at the same active geometry and budget. If suitable QM reference data
are available, separate prescribed-field response from accuracy for a real coordinated ion. Do not
make a new QM dependency or large protein benchmark project a requirement for this software change;
without reference data, state that quantitative scientific accuracy remains unvalidated.

Review the entire flow for default-off compatibility, active-only coverage, fixed values, charge
accounting, target isolation, and provenance. Update owning docs only with implemented behavior, and
leave unfinished scientific qualification in TODO.md. Record final commits and validation summary here.

## Progress Notes

- 2026-10-03: Read the research handoff and inspected integration points. Created `fragments` and
  committed METALS.md as `4324b11`. Implementation and contract approval remain pending.
- 2026-10-03: Step 1 proposal prepared by Luna and reviewed by Astra; addressed review fixes in
  `PLAN.md`. User contract approval remains pending; no implementation step is complete.
- 2026-10-03: User clarified model responsibilities: Luna only for coding; Astra for architecture,
  critical decisions, and review. Astra rechecked the current proposal, clarified the public method-input
  boundary and source lifetime/validation requirements, and retained the existing budget/partition
  design. The earlier Luna design pass is historical, not the workflow for subsequent steps. User
  contract approval remains pending.
- 2026-10-03: User authorized starting with EEM but requested smaller reviewable boundaries, such as API
  preparation. Astra scoped slice 2a to the native request types and fail-closed assessment handling;
  subsequent API and execution work remains outside this slice.
- 2026-10-03: Slice 2a implemented and reviewed. Absent/empty embedding preserves ordinary assessment;
  nonempty sources throw before planning. Stopped for user review without starting partitioning or
  numerical execution. Step 2 remained incomplete at that point.
- 2026-10-03: Generalized the API name to `FixedChargeEmbedding` / `fixed_charge_embedding` at user
  request. Indexed atomic sources can also represent future frozen components; initial isolated-source
  scope and current fail-closed behavior are unchanged. Fragment validation and component-total metadata
  remain optional future work. GCC debug full suite (56/56) and focused Clang planning test passed.
