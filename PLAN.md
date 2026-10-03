# Fixed-Charge Embedding Implementation Plan

Status: EEM implementation authorized, beginning with smaller independently reviewable slices.
The request, validation, partition, method-input, EEM numerical, and assessment slices (2a-2g) are
committed. Private result reconstruction (2h) is implemented and reviewed; facade embedding execution remains disabled.
Branch: `fragments`. Research baseline: `4324b11` ([METALS.md](METALS.md)).

## Goal and Scope

Calculate active-molecule charges in the electrostatic field of explicitly prescribed, fixed in-target
ion charges, without requiring fitted parameters for those ions. Removing ions and appending their
charges after an unperturbed calculation does not meet this goal.

The initial deliverable covers EEM, SQE/SQE+q0/SQE+qp, and QEq, followed by supported cutoff/cover
execution and matching native, Python, and CLI behavior. Start with isolated monatomic sources. Keep
parameter-covered ligands active. ABEEM and EQeq are bounded follow-ons; EQeq+C and frozen molecular
fragments are optional, not release blockers.

No new fitted parameter sets, automatic protonation, covalent cutting/capping, fragment-charge generator,
preset library, independent external charge clouds, or response variants for every method.

This file is the user-requested implementation progress tracker. Keep research rationale in METALS.md
and implemented contracts in the owning `docs/` files. Do not duplicate this checklist in TODO.md.

## Working Rules

- **GPT-6 Astra owns architecture, public API design, scientific/coupling decisions, and code review.**
  Astra resolves critical choices before delegating a coding slice; an agent review does not substitute
  for user approval where this plan requires it.
- **GPT-6 Luna is used only for coding**, including focused tests and implementation-driven documentation
  updates within Astra's agreed design. If coding exposes an architectural or scientific ambiguity,
  Luna reports it to Astra rather than choosing a new contract or expanding scope.
- For each slice, Astra defines the boundary, Luna implements and runs the relevant checks, and Astra
  reviews correctness, scientific conventions, ownership/mapping, and regressions. Luna addresses review
  findings; Astra verifies the fixes before reporting the slice ready for user review. Stop at that
  boundary instead of starting the next step automatically.
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
  preparation. The EEM request API slice below is authorized; QEq coupling remains pending approval.

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

Use the original supplied formal-charge sum as `Q_original`; for each affected target,
`Q_active = Q_original - sum(fixed charges)`. Imported `core::Atom::formal_charge()` is an integer
defaulting to zero with no presence flag, so a fixed Mg charge of +2 does not repair an imported formal
charge of zero. Callers must intentionally prepare the target's formal charges when they need a
nonzero original total. Do not infer charge from element/component ID, apply CCD rules, or change the
target total based on prescribed source charges.

### Contract Decisions Proposed For Approval

| Decision | Proposed contract |
| --- | --- |
| Scope and selection | Sources belong to the same input molecule/target as active atoms; no cross-record or external point cloud. Each source is an explicitly indexed atom with a finite real prescribed charge. Initial support is isolated atoms only: reject duplicate indices, out-of-range indices, selected atoms with any graph bond, and a target with no active atoms. An empty source list normalizes to embedding disabled before capability checks and produces no embedding provenance. No element-based inference. |
| Coordinates and identity | Use original target conformer coordinates for each calculation target. A source index identifies the same atom across conformers; evaluate its position in the corresponding conformer. Preserve original molecule, conformer, and atom order in returned results. |
| Budget | Preserve the original supplied formal-charge sum as `Q_original`; calculate `Q_active = Q_original - sum(fixed charges)`. Fixed charges are included exactly once in reconstructed output. A prescribed Mg +2 does not repair an imported formal-charge total of zero; callers must prepare formal charges intentionally. SQE and SQE+q0 reject budgets incompatible with their existing component/seed constraints; never renormalize them. |
| Validation | Reject nonfinite charges/totals/coordinates, coincident active/source sites, and invalid graph/source scope. Do not clamp distances. Do not define a universal near-contact radius or warning absent a method-specific scientific basis. |
| Ownership | `AssessmentRequest` owns raw selectors. `AssessmentResult` owns its validated partition, source values, mappings, and geometry lifetime alongside its current molecule/prepared-feature owners. Its facade execution passes internally validated target-local context tied to the active prepared data and candidate. `CalculationInput` receives a read-only non-owning `std::span<const FixedPointSource>` (each source has `core::Position position` and `double charge`), valid throughout the method call. Coupling stays method-specific. |
| Lower-level execution | Keep the existing public `CalculationRequest` unchanged and non-embedding initially: it consumes already prepared/classified data and cannot accept raw selectors or repartition. Only the owned assessment facade partitions, validates, executes with its tied internal context, and reassembles results at the facade boundary. |
| Capability/resource | Add an explicit embedding capability to method requirements/assessment, independent of `full`/`cutoff`/`cover`; initially only the methods implemented in steps 2-4 advertise it, and only full execution until step 5. Unsupported candidates/modes are rejected, never silently downgraded. Resource assessment counts active solve size plus `N_active * N_source` source-field work using the existing complexity/resource mechanism; no new planner. |
| Provenance | Add structured embedding provenance to effective/result output: source molecule/atom indices and exact prescribed charges, caller-supplied charge-provenance label, source kernel convention, per-target original and active totals. Keep it separate from execution mode/radius. The existing native `ExecutionResult` and adapter `ChargeCalculationResult` are the result boundaries; JSON serializes the structure, molecular charge arrays remain original-order. |
| QEq | Proposed, not settled: analytic infinite-source-hardness limits from METALS.md 4.6, separately for each overlap option and without passing infinity to the finite-parameter kernel. Keep active-active behavior unchanged and identify this convention in provenance. This scientific convention requires explicit approval before step 4. |

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
charge, and provenance. Both use original supplied formal-charge sums. Exact Python/CLI spelling remains
for their implementation steps and is not a step-1 selector-syntax approval; neither may parse chemistry
or create its own partition policy.

Review outcome: user authorized starting EEM implementation, with API preparation as a smaller review
boundary. Implement the contract incrementally without enabling incomplete execution paths. QEq coupling
still requires explicit approval before step 4. For SQE and SQE+q0, incompatible active budgets reject
rather than changing seed policy.

## 2. Deliver EEM Full Execution End to End

- [ ] Implement one complete native path from assessment through reconstructed output.

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
- [x] **2c: Geometry and charge-budget validation.** For each affected molecule, require conformers and
  finite coordinates for every atom in every conformer; reject exact active/source coincidence without
  a distance clamp or near-contact threshold. Validate finite prescribed-source sums and the budget
  `Q_active = Q_original - sum(fixed charges)`, using the original supplied formal-charge total. Keep
  source selection, geometry, and totals target-local; do not infer formal charges or require prescribed
  charges to equal them. Valid nonempty requests still reach the unsupported-execution gate. No owned
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
  metadata preservation and verified full GCC regressions. Commit: pending user review.
- [ ] **2i and later:** define the next small boundary after review of 2h, working toward the end-to-end
  requirements below. Keep each intermediate state fail-closed.

Primary files: `include/chargefw/calculation/{assessment,calculation}.h`,
`src/calculation/{assessment,calculation,full_execution,target_execution}.*`,
`include/chargefw/methods/{calculation_input,method_requirements}.h`, and
`src/methods/builtin/eem.cpp`.

- Prepare the partition once before candidate checks. Own active molecules and their features with safe
  lifetimes; retain disconnected active components together and preserve original conformer identity.
- Declare embedding capability explicitly. Initially enable EEM/full only; automatic planning must not
  choose an unsupported method or reduced mode, or override resource policy to manufacture a plan.
- Enforce embedding compatibility in candidate assessment before classification and in existing
  execution validation before dispatch. Direct method callers are responsible for declared preconditions;
  do not add guards to algorithms that do not implement embedding or introduce a dispatch wrapper.
- Assess prerequisites and classify active atoms/bonds only with the unchanged matcher. Translate active
  diagnostic indices back to original atom/bond indices.
- Pass the active budget to the solve and subtract `sum(kappa * fixed_charge / distance)` from its RHS,
  using the selected EEM set's `kappa`. Do not read ion atom parameters.
- Reassemble original-order charges once at the common result boundary. Fixed values are copied exactly;
  reusable plans, progress, cancellation, and result identity retain their existing guarantees.
- Add structured embedding provenance and account for active solve size plus source-field work in
  execution resource assessment. A simple cost term is enough; no resource-planner rewrite.
- Document the supported native/full slice in `docs/PROJECT.md` and `docs/NATIVE.md`.

Focused check: one small active pair plus a missing-parameter ion, compared with an independently
assembled constrained EEM system. Demonstrate nonzero active response, exact fixed charge, conserved
total, and reusable original-order output. Add a compact invalid-request table to existing planning
tests for crossing bonds, invalid selections/geometry, unsupported modes, and missing active coverage.

## 3. Add the SQE Family

- [ ] Enable full embedding for SQE, SQE+q0, and SQE+qp through the shared solver.

Primary files: `src/methods/builtin/{sqe,sqeq0,sqeqp}.cpp` and method prerequisite checks.

- Subtract the external potential before incidence projection. Preserve the existing `-H*s + D*s`
  seed correction and all active-active/bond-hardness terms.
- Couple each point source using the active Gaussian width:
  `erf(r / (sqrt(2) * abs(width))) / r`, with `1/r` for zero width. No extra Coulomb factor.
- Preserve ordinary SQE component neutrality and SQE+q0 formal seed totals. Reject incompatible budgets.
- Normalize SQE+qp fitted seeds once over the whole active target to `Q_active`; fixed sources neither
  require fitted `q0` nor participate in normalization. Do not switch to per-component normalization.
- Update method capability and the shared/native documentation for this supported slice.

Focused check: extend existing `test_sqe`/`test_sqeqp` fixtures with a prescribed source and a small
explicit RHS reference, including charged SQE+qp normalization and SQE+q0 budget rejection. Parameterize
the kernel edge checks rather than creating a separate test suite for each family member.

## 4. Add QEq With the Approved Source Kernel

- [ ] Enable QEq/full after the step 1 coupling decision is recorded.

Primary file: `src/methods/builtin/qeq.cpp`.

Keep active-active overlap behavior unchanged. Implement only the approved source coupling, subtract
its potential, and use the active budget. If adopting the proposed limits, preserve Nishimoto-Mataga-
Weiss's `17.28/r`, Ohno-Klopman's active-hardness softening, and `14.4/r` for the other four options.
Validate the active-hardness domain required by the limit and identify the convention in provenance.

Focused check: one table-driven kernel/reference-solve test across the supported overlap options,
reusing the EEM fixture and budget/mapping infrastructure. Update shared/native method documentation.

## 5. Propagate Embedding Through Cutoff and Cover

- [ ] Enable reduced embedding only where its method-specific path is complete.

Primary files: `src/calculation/{reduced_execution,cutoff_execution,cover_execution}.cpp`.

- Build references and conserved groups on the active graph with the active target budget, not its
  unadjusted formal total. Preserve the one-time SQE+qp normalization from step 3.
- Supply the complete fixed-source environment to every fragment, including sources beyond its radius.
  Start by reusing direct source evaluation; precompute/project potentials only if necessary, without
  introducing a cache framework or changing coupling conventions.
- Apply conservation corrections to active results before source reinsertion. Never include fixed
  sources as fragment variables, correct their values, or double-count their contribution.
- Preserve source-index diagnostics and finite-radius cover's existing order-dependent pivot behavior.
- Update execution capability/resource handling and `docs/PROJECT.md` for supported combinations.

Focused check: extend `test_reduced_execution` with a distant source outside the fragment radius and
whole-active-radius agreement with full execution. Exercise cutoff and cover with EEM and SQE+qp to
cover ordinary budgets and normalized component references; reuse method-level kernel tests.

## 6. Expose the Policy in Python

- [ ] Support the same request and provenance through assessment, direct calculation, and reusable plans.

Primary files: `python/chargefw/{calculation,_calculation_options}.py`, `python/src/calculation.cpp`,
associated value types, exports, and extension stubs.

Translate the approved molecule/atom-index request into native data. Keep immutable request snapshots,
plan/collection checks, conformer handling, and adapter source mappings. Python must not partition the
graph or implement its own charge/capability policy. Update `docs/PYTHON.md`.

Focused check: one public workflow exercising assessment, direct calculation, and plan reuse over a
multi-conformer input, with matching native results and provenance. Run existing calculation/output
and typing checks rather than duplicating numerical tests in Python.

## 7. Expose the Policy in the CLI and Serialized Results

- [ ] Add explicit source input to both `calculate` and `applicability`, with matching output provenance.

Primary files: `apps/chargefw/{input,output}.cpp`, shared result serializers and their schemas.

Use unambiguous molecule/atom scope and index base for the selector syntax settled during CLI
implementation. Parsing translates into the same native assessment request; do not add a component-
selection language or infer Mg2+ from PDB/mmCIF element names. Retain normal charge arrays and original
atom mappings in every output path. Update `docs/CLI.md`, `docs/FORMATS.md`, and one concise executable
usage example.

Focused check: one installed-CLI round trip for applicability and calculation, plus malformed-input
coverage in existing argument tests. Check structured provenance against the output schema. Use CTest
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

## Optional Follow-ons

These are separately reviewed extensions, not prerequisites for completing steps 1-8.

- [ ] **ABEEM:** drive both atomic and covalent-radius-weighted bond-center sites with the set's `k/r`;
  validate source distances to both site types and retain half-bond output redistribution. One small
  atom/bond-site reference solve is sufficient to exercise the new coupling.
- [ ] **EQeq:** approve finite elemental screening versus the parameter-free `8.64/r` limit before
  implementation; reuse the existing solve and embedding checks.
- [ ] **Explicit molecular sources:** allow complete mapped per-atom vectors for wholly frozen graph
  components, with declared total/provenance and crossing-bond rejection. Internal frozen bonds do not
  enter the active solve. Add one molecular-vector example; do not add charge generation or presets.
- [ ] **EQeq+C, only if requested:** approve a correction-boundary variant first. Active-active-only
  correction is the simplest candidate; overwriting final fixed values after ordinary correction is
  not charge-conserving and is not acceptable.

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
  subsequent API and execution work remains outside this slice. QEq coupling approval remains pending.
- 2026-10-03: Slice 2a implemented and reviewed. Absent/empty embedding preserves ordinary assessment;
  nonempty sources throw before planning. Stopped for user review without starting partitioning or
  numerical execution. Step 2 remains incomplete.
- 2026-10-03: Generalized the API name to `FixedChargeEmbedding` / `fixed_charge_embedding` at user
  request. Indexed atomic sources can also represent future frozen components; initial isolated-source
  scope and current fail-closed behavior are unchanged. Fragment validation and component-total metadata
  remain optional future work. GCC debug full suite (56/56) and focused Clang planning test passed.
