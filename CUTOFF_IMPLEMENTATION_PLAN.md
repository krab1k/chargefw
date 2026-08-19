# Cutoff, Applicability, and Execution Planning

## Status and purpose

This document is a review plan, not a frozen specification. It records the current direction for
introducing cutoff and later cover execution without coupling those policies to individual method
options or duplicating selection logic in the CLI. Each implementation step is intended to be a
small, independently reviewable commit. Details may be revised after a commit is reviewed or after
compatibility and benchmark results provide better evidence.

The first objective is a correct, explicit, serial reference implementation. Efficient spatial
search, cover execution, and parallel scheduling come later. Full calculations must remain unchanged
unless a caller explicitly selects another execution mode.

This plan complements, but does not replace, the ownership of `PROJECT.md`, `TODO.md`, `README.md`,
and `AGENTS.md`. Implementation commits must update the owning document whenever they change current
architecture, completed work, CLI usage, or tracked deliverables.

## Agreed direction

1. Keep scientific method applicability, execution availability, selection, and calculation as
   distinct concepts, while retaining a facade that composes them for ordinary callers.
2. Let applicability discovery report that a method/parameter-set pair can be scientifically valid
   even when full execution is discouraged by the shared size threshold.
3. Represent `auto` only as a caller selection preference. Every selected calculation and every
   result must contain a concrete effective mode: `full`, `cutoff`, or `cover`.
4. Use one application-level atom-count threshold for expensive full calculations. Do not put
   method-specific atom thresholds in `MethodRequirements`.
5. Treat that threshold as an automatic-selection safeguard, not a hard scientific limitation.
   Explicit `full` overrides it and produces a warning before solver allocation.
6. Allow the shared threshold to be configured, including an unlimited value. With automatic
   execution and an unlimited threshold, full remains size-eligible for every otherwise valid
   candidate.
7. Require finite cutoff and cover radii of at least 8 angstrom. The eventual automatic/default
   radius must be explicit in provenance and may be revised after compatibility and benchmark work.
8. Never silently change execution mode, radius, classification policy, parameter set, method, or
   charge correction. Automatic choices must be returned and serialized.
9. Apply permissive parameter typing during applicability only. Classification must happen once; a
   selected candidate carries the resolved classifications into planning and calculation without
   carrying or reinterpreting the original classification options.
10. Start with deterministic serial code and a simple neighbor scan. Do not add a spatial-index
    dependency, parallelism, or cover abstractions before the reference cutoff behavior is tested.
11. Use `cutoffnew` as research material only. Do not merge or cherry-pick it wholesale.

## Adopted public workflows

The advanced/native workflow separates policy resolution from execution:

```text
ApplicabilityRequest
  -> ApplicabilityResult
  -> caller selection or select_applicable_method()
  -> CalculationRequest { prepared molecules, selected ApplicableMethod }
  -> CalculationResult
```

`ApplicabilityRequest` owns classification policy. `ApplicableMethod` owns the resulting atom/bond
classifications. Execution-only `CalculationRequest` consumes the selected candidate and does not
accept classification options.

The convenience workflow composes the same stages:

```text
ApplicationCalculationRequest
  -> preparation
  -> applicability
  -> deterministic selection
  -> execution
  -> ApplicationCalculationResult
```

The application request may expose classification, selection, resource, and execution preferences
because the facade resolves them internally. It must delegate to the same applicability, selection,
and execution operations rather than duplicate their behavior. As execution planning is introduced,
a concrete mode or internal immutable plan will be inserted between selection and execution; an
unresolved `automatic` preference will never be passed directly to calculation.

## Non-goals for the initial cutoff work

- No silent ChargeFW2-style atom-count switching.
- No automatic chemistry preparation or topology repair.
- No per-method atom thresholds.
- No hardware-specific memory estimator in the first implementation.
- No promise that an explicitly forced full calculation can recover from operating-system OOM
  termination.
- No cutoff or cover claim for SQE, SQE+q0, or SQE+qp until their fragment semantics are specified and
  validated.
- No optimized KD-tree, fragment cache, task scheduler, OpenMP integration, or nested solver thread
  management in the reference implementation.
- No broad CLI rewrite unrelated to calculation selection and provenance.

## Terminology and conceptual model

### Scientific applicability

A method/parameter-set pair is scientifically applicable when the molecules satisfy requirements
that cannot be bypassed by an execution override, including:

- required coordinates, topology, bond orders, formal charges, and element properties;
- method-specific molecular restrictions;
- parameter-set association and required parameter names;
- successful strict or permissive parameter classification, according to the request.

Missing scientific prerequisites remain hard failures. `--execution full`, an unlimited threshold,
or any future unsafe option must not bypass them.

### Execution availability

For a scientifically applicable pair, execution assessment reports the status of each mode:

- `full` is implemented for all current methods, but automatic use may be discouraged above the
  shared threshold for methods declared to have expensive full complexity;
- `cutoff` is available only after a method has a tested cutoff implementation;
- `cover` is available only after a method has a tested cover implementation;
- cutoff and cover expose radius constraints rather than pretending that an unspecified radius is a
  concrete executable policy.

Exceeding the shared full threshold is not scientific inapplicability. It is an execution warning and
an automatic-selection constraint that explicit full execution may override.

### Selection preference and concrete execution policy

The request may contain an execution preference such as `automatic`, `full`, `cutoff`, or `cover`.
The selected calculation must contain a concrete execution policy. `automatic` must never survive as
the effective mode in a result.

### Calculation plan

A calculation plan is an executable, immutable choice containing at least:

- method and optional parameter set;
- validated method options;
- strict/permissive classification policy;
- classifications produced under that policy;
- concrete execution mode and optional radius;
- resource warnings and the effective shared threshold;
- enough provenance to reproduce and explain selection.

The implementation does not need to expose all of these as one large public type immediately. The
important invariant is that execution consumes the exact validated choices produced by planning and
does not independently reinterpret request defaults.

## Minimal public-policy model

The exact names may change during review. Prefer small value types over a general policy framework.
A likely initial shape is:

```cpp
inline constexpr std::size_t default_full_atom_threshold = 20'000;
inline constexpr double minimum_reduced_radius = 8.0;

enum class ExecutionMode {
    full,
    cutoff,
    cover,
};

enum class ExecutionSelectionKind {
    automatic,
    full,
    cutoff,
    cover,
};

struct ExecutionSelection {
    ExecutionSelectionKind kind = ExecutionSelectionKind::automatic;
    std::optional<double> radius;
};

struct ExecutionPolicy {
    ExecutionMode mode = ExecutionMode::full;
    std::optional<double> radius;
};

struct ResourcePolicy {
    // A value is the shared threshold. nullopt means unlimited.
    std::optional<std::size_t> full_atom_threshold = default_full_atom_threshold;
};
```

Rules:

- `full` rejects a supplied radius because the radius has no meaning for full execution;
- explicit `cutoff` and `cover` require a finite radius of at least 8 angstrom;
- automatic reduced execution uses an explicitly defined default radius only after that default is
  accepted for the reference implementation;
- `ResourcePolicy` is application/planning policy, not a method option;
- explicit `full` remains executable above the threshold and adds a warning;
- automatic selection treats expensive full execution above the threshold as discouraged;
- an unlimited threshold removes only the atom-count safeguard, not scientific prerequisites or
  solver errors.

The old implementation used a 12-angstrom default. The prototype branch used 8 angstrom. Eight is
accepted as the minimum, but the automatic default should be confirmed at the planning-selection
commit. Until then, tests should not accidentally establish an unsupported scientific recommendation.

## Shared expensive-full rule

`MethodRequirements` should retain complexity metadata and execution capabilities but no atom-count
value. A central helper in calculation planning determines whether the shared full threshold applies.
The initial rule should be deliberately simple and testable: full execution is expensive when the
declared time complexity is cubic or the declared memory complexity is quadratic in atoms, bonds, or
their combination.

The helper must explicitly enumerate the relevant `ComplexityTerm` values rather than depend on enum
ordering. If review shows that the existing complexity terms are insufficient, adjust the terms in a
separate focused change rather than adding method-specific thresholds.

For a collection-level plan, compare the threshold to each molecule independently. If any molecule
exceeds it, a collection-wide full plan receives the warning or is excluded from automatic selection.
The initial implementation uses one execution mode for the whole collection; it must not silently mix
full and cutoff assignments within one result.

## Applicability result shape

`ApplicabilityRequest` is now established for classification policy. Commit 3 extends it with resource
policy so neither concern becomes another positional argument:

```cpp
struct ApplicabilityRequest {
    const features::PreparedMoleculeCollection& molecules;
    std::span<const Method* const> methods;
    std::span<const parameters::ParameterSet> parameter_sets;
    parameters::ClassificationOptions classification_options{};
    calculation::ResourcePolicy resources{};
};
```

Namespace dependencies must be checked during implementation. If putting `ResourcePolicy` in
`calculation` causes an undesirable dependency from `methods`, place the small policy value types in
a neutral public execution-policy header rather than duplicating them.

For each scientifically applicable method/parameter pair, report a small execution assessment:

```cpp
enum class ExecutionAvailability {
    available,
    available_with_warning,
    unsupported,
};

struct ExecutionAssessment {
    ExecutionMode mode;
    ExecutionAvailability availability;
    std::optional<double> minimum_radius;
    std::vector<PrerequisiteIssue> issues;
};
```

This is illustrative, not a requirement to use exactly these enums. Avoid introducing a hierarchy or
generic diagnostic framework. Existing prerequisite issues should be reused where practical, with a
new typed issue only when needed to distinguish unsupported execution from an advisory resource
warning.

An expensive large SQE candidate before SQE cutoff exists should be representable as:

```text
scientifically applicable
full: available with resource warning; explicit override allowed
cutoff: unsupported
cover: unsupported
```

After validated SQE cutoff is eventually added, only the cutoff assessment changes.

## Permissive parameter typing

Permissive typing is exposed by `methods::ApplicabilityRequest` and by
`calculation::ApplicationCalculationRequest`, because the application facade initiates applicability
internally. The execution-only `calculation::CalculationRequest` does not accept classification
policy. This preserves one classification decision without redesigning classification.

Required invariants:

1. Applicability uses the request's `ClassificationOptions` for every method/parameter pair.
2. The successful classifications are stored with each `ApplicableMethod` candidate.
3. Calculation constructs `ParameterView` from those stored classifications and never reclassifies.
4. Explicit method/parameter selection and automatic selection use the same classification policy.
5. CLI, future bindings, and direct native callers receive the same result for the same policy.
6. Result provenance eventually records whether permissive typing was enabled without requiring the
   policy to be retained on `ApplicableMethod` or passed back into execution.
7. A strict failure that becomes successful under permissive typing has a focused test at both the
   applicability and facade-calculation levels.

The initial change need only report whether permissive matching was enabled. Recording exactly which
atoms or bonds required fallback is useful future diagnostics but is not required to begin cutoff and
must not expand this change unnecessarily.

## Automatic and explicit selection rules

The initial selector should preserve the existing method-priority and parameter-priority ordering.
For each ranked scientifically applicable pair:

1. Explicit `full` selects full even above the shared threshold and adds a warning.
2. Explicit `cutoff(radius)` or `cover(radius)` requires method support and a valid radius; it never
   falls back to another mode.
3. Automatic selection prefers full when it is not discouraged by the shared threshold.
4. For an expensive candidate above a finite threshold, automatic selection may use cutoff and then
   cover only when the corresponding implementation exists and an accepted automatic radius is
   available.
5. If no reduced mode can be selected, continue to the next ranked scientifically applicable pair.
6. If no concrete plan can be made, return applicability and execution diagnostics rather than
   silently forcing expensive full execution.
7. With an unlimited threshold, full is never discouraged by size, so the automatic path selects
   full for the first otherwise valid ranked candidate.

The ordering of cutoff versus cover is initially cutoff before cover. There is no second hidden
threshold analogous to ChargeFW2's 80,000-atom switch. Benchmark evidence may later justify a
different explicit automatic policy.

Before automatic reduced execution is enabled in the CLI, review and settle the automatic radius.
If that decision is deferred, automatic selection should report the reduced alternatives and require
an explicit mode/radius rather than inventing a value.

## CLI direction

Add CLI options only after the application facade implements the same behavior:

```text
--method <id>
--parameter-set <id>
--execution auto|full|cutoff|cover
--radius <angstrom>
--full-atom-threshold <count|unlimited>
--permissive-types
```

Expected behavior:

- `--execution full` is the deliberate resource override; no generic `--unsafe` flag is initially
  necessary;
- if explicit full exceeds the finite threshold, print a warning to stderr before fragment or solver
  allocation and retain the warning in result provenance;
- `--full-atom-threshold unlimited` disables only automatic size discouragement and is itself
  reported;
- radius is rejected for full and required for explicit cutoff/cover;
- unsupported mode selections fail clearly and do not fall back;
- method and parameter options restrict facade candidates rather than implementing selection in the
  CLI;
- permissive typing is passed through the facade rather than used by an independent CLI classifier;
- automatic choices and effective policy are included in JSON and molecular-output metadata where
  those formats support provenance.

## Use of the `cutoffnew` prototype

Reusable ideas:

- an execution-policy header with an 8-angstrom minimum;
- fragment data in `features`, outside `core::Molecule`;
- explicit local-to-source atom and bond maps;
- induced source bonds and conformer preservation;
- projection of immutable whole-molecule parameter classification into a fragment;
- a serial one-fragment-per-center reference algorithm;
- proportional fragment target charge and final uniform correction for ChargeFW2 EEM-family parity.

Problems to correct rather than copy:

- execution policy was not considered during applicability;
- the existing large-molecule prerequisite still rejected cutoff-capable candidates;
- the stashed CLI manually filtered the method registry, duplicating facade policy;
- unsupported execution could be discovered only after candidate selection;
- any non-full policy was dispatched as cutoff, which is unsafe once cover exists;
- an internal single-target helper became public solely to connect implementation files;
- solver failures were wrapped as invalid arguments;
- policy, radius, correction, and warnings were absent from result provenance;
- there was no above-threshold test proving that cutoff changes execution availability;
- only early EEM/QEq experiments existed and no SQE semantics were defined.

## Commit-by-commit implementation plan

Every commit below must build and pass its focused tests. Do not combine later commits merely because
the code is already understood. Review may change subsequent steps.

### Commit 1: Carry permissive classification through applicability

Status: implemented. Review established the explicit and convenience-facade workflows described
below; this status does not complete later CLI or provenance deliverables.

Scope:

- introduce `ApplicabilityRequest` with the current molecules, method candidates, parameter sets, and
  `ClassificationOptions`;
- update `find_applicable_methods()` to consume the request;
- pass classification options into `check_parameter_prerequisites()`;
- retain successful atom/bond classifications with each `ApplicableMethod` candidate;
- add classification options to `ApplicationCalculationRequest`, because the convenience facade
  performs applicability internally;
- make `CalculationRequest` execution-only over a prepared collection and an already-selected
  `ApplicableMethod`; do not add classification options to it or reclassify during execution;
- expose deterministic `select_applicable_method(ApplicabilityResult)` for callers that want the
  standard priority policy without using the owned facade;
- update internal call sites without changing default strict behavior;
- do not add execution modes in this commit.

Focused tests:

- an exact match remains preferred when permissive typing is enabled;
- a strict candidate is rejected when only the permissive fallback matches;
- the same candidate is applicable with permissive typing;
- automatic facade calculation uses the stored permissive classification successfully;
- the explicit applicability/select/calculate workflow uses the same stored classification;
- default requests remain strict.

Review checkpoint:

- confirm that no calculation path reclassifies;
- confirm that classification policy is absent from execution-only `CalculationRequest` and
  `ApplicableMethod` while the resolved classifications remain available;
- confirm that the request shape is suitable for later resource policy without becoming a general
  options bag.

Likely documentation:

- update `PROJECT.md` only if the public workflow description becomes stale;
- keep the corresponding TODO item open until CLI/provenance support is complete.

### Commit 2: Add minimal execution and shared resource policy types

Status: implemented. These types validate and carry future execution preferences only; the current
facade still performs the existing full calculation until execution availability and concrete plan
selection are implemented in Commits 3 and 4.

Scope:

- add concrete execution mode/policy value types;
- add caller execution selection with `automatic` distinct from concrete policy;
- add `ResourcePolicy` with one shared default full atom threshold and unlimited representation;
- add finite/minimum-radius validation;
- add the caller execution preference to `ApplicationCalculationRequest` without changing
  calculation behavior yet;
- keep `ResourcePolicy` ready for `ApplicabilityRequest` in Commit 3, where execution availability is
  assessed;
- do not add an unresolved execution preference to execution-only `CalculationRequest`; a concrete
  selected mode or internal plan reaches execution only when Commit 4 introduces plan selection;
- document that explicit full is the threshold override;
- do not add fragment code or advertise cutoff support.

Focused tests:

- default values;
- full with no radius;
- rejection of a radius attached to full;
- rejection of missing, NaN, infinite, and sub-8-angstrom reduced radii;
- acceptance of exactly 8 angstrom and larger finite values;
- finite shared threshold and unlimited threshold representation.

Review checkpoint:

- choose the automatic reduced radius or explicitly defer automatic reduced execution;
- confirm CLI spelling before it is exposed.

### Commit 3: Make applicability report execution availability

Scope:

- remove per-method atom threshold and hard large-molecule rejection fields;
- retain method complexity and cutoff/cover capabilities;
- add the central expensive-full complexity helper;
- add `ResourcePolicy` to `ApplicabilityRequest`;
- report full, cutoff, and cover assessments for scientifically applicable candidates;
- report threshold exceedance as an advisory/automatic-selection issue, not scientific rejection;
- keep missing inputs and parameter failures as hard rejection;
- do not mark built-in methods cutoff-capable until their executor and tests exist.

Focused tests:

- inexpensive full methods remain automatically available above the threshold;
- an expensive method is full-available normally below the threshold;
- above the threshold, expensive full is reported with an overridable warning;
- unlimited threshold removes the size warning;
- a lower custom threshold permits testing the same behavior on a small fixture;
- unsupported cutoff and cover are reported without making the scientific pair disappear;
- collection behavior checks every molecule;
- strict and permissive parameter coverage produce the same execution assessments after successful
  scientific applicability.

Review checkpoint:

- inspect all built-in complexity declarations before relying on them for shared-threshold behavior;
- confirm whether the issue belongs in existing prerequisite diagnostics or a small execution-specific
  diagnostic type.

### Commit 4: Add concrete plan selection without reduced calculation

Scope:

- extend the existing deterministic candidate selector with concrete plan selection;
- create or assemble a concrete full calculation plan from an assessed candidate;
- implement explicit full override and warning behavior;
- implement explicit rejection of unsupported cutoff/cover selections;
- make automatic selection obey the shared threshold and unlimited setting;
- retain diagnostics when no concrete plan exists;
- route the existing facade through assess, select, and execute while full calculation remains the
  only implemented concrete mode.

Focused tests:

- existing method and parameter priority/tie behavior is unchanged for ordinary full calculations;
- explicit full above threshold calculates and returns a warning;
- automatic full above threshold is not silently forced for an expensive method;
- unlimited threshold makes automatic full size-eligible;
- explicit unsupported cutoff/cover fails before calculation;
- unavailable explicit method or parameter IDs still fail without fallback;
- the effective classification policy is recorded in result provenance, while execution continues to
  reuse the candidate's stored classifications;
- the selected concrete execution policy survives into the result.

Review checkpoint:

- confirm the automatic behavior when a higher-priority expensive candidate is discouraged but a
  lower-priority inexpensive full candidate exists;
- avoid stabilizing a public `CalculationPlan` type unless an external caller already needs it.

### Commit 5: Add deterministic spatial fragment support

Scope:

- adapt the prototype's `SpatialFragment` into `features`;
- preserve selected source atom order deterministically while identifying the center explicitly;
- preserve atomic number, formal charge, atom name, one source conformer, and induced source bonds;
- retain local-to-source and source-to-local atom mapping plus local-to-source bond mapping;
- project stored atom and bond classifications without reclassifying fragment environments;
- use the existing linear conformer-neighbor operation initially;
- keep single-target execution plumbing internal to the library target where possible.

Focused tests:

- radius boundary inclusion;
- center inclusion and mapping;
- atom, bond, order, formal-charge, name, coordinate, and conformer preservation;
- induced bond omission at fragment boundaries;
- atom-only, bond-only, and atom-plus-bond classification projection;
- empty and single-atom molecules where valid;
- invalid atom/conformer indices and non-finite/invalid radii;
- deterministic mapping independent of neighbor discovery order.

Review checkpoint:

- compare fragment ordering with ChargeFW2 before numerical compatibility fixtures depend on it;
- confirm that classification projection, rather than fragment reclassification, is correct for
  compatibility and induced-topology semantics.

### Commit 6: Implement serial EEM cutoff as the first executable reduced mode

Scope:

- add a central serial cutoff executor outside individual `Method` state;
- support one fragment per source atom and extract the center charge;
- reproduce ChargeFW2 EEM-family target-charge allocation for the first EEM implementation;
- apply and report the final uniform total-charge correction;
- preserve molecule, conformer, atom, method, and parameter-set identity;
- mark only EEM cutoff-capable after its focused tests pass;
- let applicability now report EEM cutoff as available;
- let explicit cutoff produce a concrete plan and execute it;
- do not implement cover, parallelism, or a spatial index.

Focused tests:

- cutoff below 8 angstrom is rejected during planning;
- EEM cutoff is available with complete parameters and unavailable otherwise;
- a custom low threshold demonstrates that cutoff can make an expensive EEM candidate automatically
  executable without constructing a 20,001-atom test molecule;
- one assignment per conformer and preserved target identity;
- atom-order preservation and finite result validation;
- total charge after correction;
- cutoff equals full when the radius contains the whole molecule, within a stated tolerance;
- deterministic repeated results;
- fragment solver failures retain molecule, conformer, and center context but remain calculation
  failures rather than invalid-request errors;
- focused ChargeFW2 compatibility fixture using identical molecule, coordinates, parameters, radius,
  target-charge policy, and correction.

Review checkpoint:

- inspect numerical deviation from ChargeFW2 before enabling additional methods;
- settle and record the automatic radius before allowing automatic cutoff in user-facing CLI paths.

### Commit 7: Expose selection, threshold, permissive types, and cutoff in the CLI

Scope:

- route method ID, parameter-set ID, classification policy, execution preference, radius, and shared
  threshold through `ApplicationCalculationRequest`;
- add the CLI options listed above;
- print resource warnings before calculation;
- fail explicit unsupported selections without fallback;
- do not filter the registry manually in `main.cpp`;
- serialize effective concrete execution mode, radius, shared threshold/unlimited state,
  threshold-exceeded warning, correction policy, and permissive-types setting;
- retain selected method and parameter IDs as today.

Focused tests:

- CLI explicit full below and above a custom low threshold;
- CLI unlimited threshold;
- CLI explicit EEM cutoff at 8 angstrom and a larger radius;
- missing/invalid radius and radius supplied with full;
- unsupported method/mode combination;
- strict classification failure and permissive success on the same fixture;
- explicit method/parameter selection and automatic selection;
- JSON provenance contains the effective concrete policy and warnings;
- exit status distinguishes invalid requests, no applicable plan, and calculation failure where the
  existing result contract permits it.

Documentation:

- update `README.md` with exact CLI syntax and defaults;
- update `PROJECT.md` current CLI and execution-policy state;
- update `TODO.md` only for deliverables fully completed by this point.

### Commit 8 and later: Expand validated EEM-family cutoff support incrementally

Use separate reviewable commits, each with method-specific full/cutoff and ChargeFW2 comparison tests:

1. QEq;
2. EQEq;
3. EQEq+C;
4. SFKEEM;
5. SMP/QEq;
6. any additional method proven to share the archived `EEMethod` behavior.

Do not mark a method cutoff-capable merely because it can technically run on a fragment. Confirm its
target-charge handling, required parameters, correction, numerical behavior, and archived support.

### Later milestone: Cover reference implementation

After cutoff compatibility is stable:

- reproduce archived pivot selection and overlap semantics deterministically;
- keep the same shared full threshold; do not introduce a hidden 80,000-atom mode switch;
- report cover radius, pivots/coverage diagnostics, overlap reconciliation, and final correction;
- validate every atom receives at least one assignment;
- add explicit and automatic plan tests only after cover behavior is numerically validated;
- optimize and parallelize only after the serial reference passes compatibility tests.

### Later milestone: SQE-family reduced execution

Treat SQE, SQE+q0, and SQE+qp as a new approximation design, not ChargeFW2 parity. Before code, write
and review the semantics for:

- cut bonds and induced fragment transfer variables;
- connected components at fragment boundaries;
- fragment target charge;
- SQE zero-component-charge behavior;
- SQE+q0 formal initial charges;
- SQE+qp parameterized initial charges and their correction;
- overlap or center reconciliation;
- final total-charge correction;
- convergence toward full as radius increases.

Once accepted, the same applicability, permissive-classification, planning, threshold, CLI, and
provenance path should support SQE without a facade redesign.

## Validation strategy

### Narrow checks first

For every commit:

1. build the directly affected target;
2. run the narrowest corresponding CTest target;
3. run `ctest --preset gcc-debug` after focused tests pass;
4. use the sanitizer preset for fragment ownership/mapping changes and before considering the
   reference implementation complete;
5. avoid unrelated formatting and inspect the diff before each commit.

### Core invariants

Every execution mode must preserve or validate:

- output charge count equals source atom count;
- all charges are finite;
- assignment count and conformer identity are correct;
- molecule and atom ordering are unchanged;
- method and optional parameter-set IDs are owned by the result;
- total-charge correction is explicit and tested;
- strict/permissive policy is the same at applicability and calculation time;
- effective execution mode and radius are concrete and reported;
- repeated serial execution is deterministic;
- unsupported modes fail before numerical calculation;
- explicit full threshold override warns before large allocation.

### Scientific and performance validation after the reference path

- compare full and cutoff over increasing radii and verify convergence;
- compare against ChargeFW2 for archived EEM-family modes;
- record per-method tolerances and intentional deviations;
- establish benchmark corpora before changing automatic thresholds or radius defaults;
- measure runtime and memory before choosing a spatial index or scheduling design;
- validate parallel output against the serial reference before enabling parallel execution by default.

## Risks and controls

### Public API churn

`ApplicabilityRequest`, execution assessment, and application request fields affect public headers.
Keep each type small, avoid inheritance, and review names before CLI/binding stabilization. Preserve a
thin compatibility overload temporarily only if downstream use requires it; do not maintain two
independent implementations.

### Applicability becoming too complicated

Do not create a generic rule engine. Continue using existing prerequisite and classification code,
then append three concrete execution assessments. Separate hard scientific rejection from advisory
resource warnings using the smallest typed distinction that works.

### Automatic selection changing scientific results

Default full behavior must remain explicit until automatic reduced radius and fallback are reviewed.
Whenever automatic reduced selection is enabled, serialize the reason, effective mode, radius, and
correction. Explicit selections never fall back.

### Out-of-memory behavior

Explicit full above the threshold may be killed by the operating system rather than throwing a C++
exception. Emit warnings before solver allocation and document that the override removes a safeguard;
do not promise recovery.

### Prototype behavior accidentally becoming specification

Use prototype fragment and correction code only after comparing it with the archived implementation.
Tests should establish intentional semantics, not merely preserve the prototype output.

### Scope expansion

Do not combine cutoff with spatial optimization, parallelism, cover, SQE semantics, broad result-schema
redesign, or unrelated CLI work. Add only the provenance fields required to report the policy actually
implemented in the current commit.

## Completion criteria for the initial cutoff milestone

The initial milestone is complete when:

- strict and permissive classification are both available through applicability, facade, and CLI;
- applicability reports scientific validity separately from full/cutoff/cover availability;
- one shared configurable full atom threshold governs automatic planning for expensive methods;
- explicit full overrides the threshold with a pre-allocation warning;
- unlimited threshold makes full always size-eligible;
- EEM cutoff with radius at least 8 angstrom is available through the native and application facades;
- the CLI can explicitly select method, parameter set, execution mode, radius, threshold, and
  permissive typing without duplicating policy;
- result provenance records effective execution and classification policies;
- serial EEM cutoff preserves mapping, conformers, total charge, and deterministic output;
- EEM cutoff has focused full-convergence and ChargeFW2 compatibility coverage;
- all debug tests pass and relevant sanitizer tests pass;
- `PROJECT.md`, `TODO.md`, and `README.md` accurately reflect the implemented state.

Efficient neighbor search, parallel cutoff, cover, and SQE-family reduced execution remain separate
follow-up milestones.
