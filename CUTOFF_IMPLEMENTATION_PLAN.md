# Reduced Execution Design Record

This document retains the design decisions and unfinished validation for cutoff and future cover
execution. Implemented state is summarized in [PROJECT.md](PROJECT.md); actionable work is tracked only
in [TODO.md](TODO.md). The former commit-by-commit plan was removed after its applicability, policy,
planning, fragment, executor, CLI, provenance, and SQE-cutoff stages were implemented.

## Current status

Implemented:

- scientific applicability separated from execution availability;
- strict/permissive classification carried once into `ApplicableMethod`;
- concrete `full`, `cutoff(radius)`, and `cover(radius)` policy values plus caller `auto` selection;
- one shared configurable full atom threshold (20,000 by default; unlimited supported);
- deterministic concrete-plan selection and explicit full override warnings in results;
- source-ordered `features::SpatialFragment` with induced bonds and classification projection;
- one-fragment-per-center cutoff execution through ordinary `Method::calculate()`, scheduled in
  parallel over independent targets or fragment centers without nested worker pools;
- one `SpatialFragmentBuilder` and nanoflann KD-tree per source conformer for cutoff fragment neighbor
  queries;
- automatic cutoff at 12 Å when an expensive full candidate exceeds the threshold;
- CLI selection and JSON requested/effective provenance;
- cutoff capability for ABEEM, EEM, EQeq, EQeq+C, QEq, SQE, SQE+q0, and SQE+qp.

Not implemented or incomplete:

- cover execution and overlap semantics;
- pre-allocation emission of explicit-full warnings (warnings are currently retained after execution);
- ChargeFW2 cutoff compatibility fixtures and supported numerical tolerances;
- broad multi-radius, charge-state, disconnected-system, performance, and memory validation;
- complete SQE-family accuracy envelopes;
- full provenance in non-JSON output.

## Invariants

1. Scientific applicability, execution availability, selection, and numerical execution remain
   distinct operations.
2. Missing scientific prerequisites are hard failures. A full override or unlimited threshold cannot
   bypass missing geometry/topology/formal charges/parameters or method restrictions.
3. `auto` is only a request preference. Every executable plan and result contains `full`, `cutoff`, or
   `cover`.
4. The shared atom threshold guides automatic execution; it is not a scientific limit. Explicit full
   may override it and must report a warning.
5. Classification occurs during applicability. Planning and execution reuse stored classifications and
   never classify fragments independently.
6. Execution mode, radius, correction, method, parameter set, and classification policy are explicit
   and reproducible; unsupported explicit choices never fall back.
7. One mode applies to the entire collection. Threshold assessment examines each molecule separately
   and must not silently mix full and cutoff assignments.
8. Full calculations remain faithful to their cited methods. Approximate capability is enabled only
   after method-specific tests and validation.
9. Source atom order, molecule/conformer identity, finite values, charge count, and deterministic
   serial output are preserved at every boundary.

## Policy model

Public constants and values live in `include/chargefw/calculation/execution_policy.h`:

```text
minimum reduced radius       8 Å
default automatic radius    12 Å
default full threshold  20,000 atoms
selection kinds        auto, full, cutoff, cover
concrete modes               full, cutoff, cover
correction policies          none, uniform
```

Rules:

- full rejects radius and charge correction;
- explicit cutoff/cover requires a finite radius of at least 8 Å;
- automatic may accept a valid radius override but not a correction override;
- reduced execution defaults to uniform correction;
- unlimited threshold removes only atom-count discouragement;
- explicit full above threshold remains available with resource issues;
- a full calculation is considered expensive when declared cubic in time or quadratic in memory,
  explicitly enumerating atom, bond, and atom-plus-bond complexity terms.

Automatic planning examines candidates in deterministic method/parameter priority order. For each
candidate it selects non-discouraged full, then supported cutoff, then supported cover. If no mode is
available it continues to the next candidate. It never silently forces discouraged full execution.

## Fragment and executor design

`features::SpatialFragment` is method-neutral. It contains:

- one induced molecule and source conformer;
- selected atoms in source order and an explicit center-local index;
- induced source bonds;
- local-to-source atom/bond mappings and source-to-local atom lookup;
- projected atom/bond classifications.

Cutoff ownership is outside the fragment. The executor:

1. builds one radius fragment per source atom using the source conformer's KD-tree;
2. prepares the fragment and projects the selected whole-molecule classification;
3. computes the method-declared fragment target charge;
4. calls the selected `Method` through normal `CalculationInput` without method-ID switches or
   downcasts;
5. validates and retains only the mapped center charge;
6. applies the selected final correction to source-ordered results.

The dense source-to-local mapping is convenient but redundant. Reassess it before retaining many
fragments or parallelizing; a builder-local lookup may be preferable. Do not optimize before the
serial reference behavior is covered by compatibility and convergence tests.

## Method semantics

### EEM/QEq-like methods and ABEEM

The generic supported policy assigns fragment target charge proportional to fragment/source atom count
and source formal charge, then optionally restores source formal charge with uniform final correction.
Whole-molecule-radius tests establish agreement with full execution for the eight declared methods,
but only EEM/QEq-like ChargeFW2 comparisons can establish archived parity.

SFKEEM remains unsupported. Its rapidly decaying kernel and global chemical-potential constraint make
neutral local-fragment targets unsuitable; forcing it through the generic executor produced poor
convergence. It needs a dedicated global sparse/truncated-kernel solver or a separately specified and
validated target policy.

SMP/QEq also remains unsupported until its archived behavior, target charge, parameters, and numerical
convergence are validated.

### SQE family

SQE cutoff is a new approximation, not ChargeFW2 parity:

- cut bonds are omitted; transfer variables exist only on induced fragment bonds;
- SQE uses zero target charge in every fragment and for final correction;
- SQE+q0 projects source formal initial charges and induced transfers preserve component totals;
- SQE+qp projects parameterized `q0`, corrects it to fragment formal charge, then applies induced-bond
  transfers;
- only center charges contribute to the source result;
- uniform final correction targets zero for SQE and source formal charge for SQE+q0/qp.

Small neutral and charged whole-molecule-radius fixtures reproduce full execution. Remaining work must
cover disconnected charged SQE+qp and multi-radius error/convergence behavior. SQE+q0's preliminary
large-structure deviation is high enough that no accuracy claim or automatic recommendation should be
made beyond the generic threshold policy.

## Cover design constraints

Cover should reuse `SpatialFragment`, not introduce a second molecule type. Execution-owned work items
must define pivots, a solve halo, retained interior atoms, source mapping, and deterministic overlap
reconciliation.

Before declaring any method cover-capable, tests must prove:

- deterministic pivot and contribution selection;
- every source atom receives at least one assignment;
- boundary atoms provide context without automatically contributing output;
- overlap estimates are reconciled deterministically;
- target-charge and final-correction behavior is explicit;
- full/cover convergence and charge conservation hold over a maintained corpus;
- diagnostics/provenance record radius, coverage, overlap, and correction.

There is no hidden second threshold analogous to ChargeFW2's 80,000-atom switch. Cover may be added to
automatic selection only after serial numerical validation.

## Validation sequence

For each reduced-execution change:

1. build the directly affected target;
2. run its focused CTest executable;
3. run `ctest --preset gcc-debug`;
4. run `clang-asan` for fragment ownership/mapping and before completing a reduced milestone;
5. compare source mapping, assignment cardinality, conformer identity, finite values, total charge, and
   repeated serial output;
6. add full-versus-reduced radii and ChargeFW2 comparison data where scientifically applicable;
7. inspect the diff and update `PROJECT.md`/`TODO.md` only when their owned state changes.

`tests/fixtures/corpus/cif/10aw.cif` and `1ek9.cif` are useful manual development structures, but their
one-off measurements are not benchmarks, tolerances, compatibility claims, or evidence for changing
the 12 Å automatic default. A maintained corpus must span methods, sizes, radii, charge states,
disconnected systems, and conformers before supported error envelopes or policy changes are claimed.
