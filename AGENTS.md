# ChargeFW Development Guide

## Required reading and document ownership

Read this file before changing the repository. Then read only the relevant sections of:

- [PROJECT.md](PROJECT.md): implemented architecture, capabilities, compatibility state, and product
  direction;
- [TODO.md](TODO.md): unfinished deliverables and their acceptance criteria;
- [README.md](README.md): supported build, test, installation, and CLI usage;
- [CUTOFF_IMPLEMENTATION_PLAN.md](CUTOFF_IMPLEMENTATION_PLAN.md): reduced-execution design decisions
  and remaining cutoff/cover validation work.

Keep ownership clear. Update the owning document when a change makes it stale; do not copy the same
status or instructions into several files. Completed work belongs in `PROJECT.md`, not as checked
history in `TODO.md`. `README.md` should remain a concise user entry point.

## Repository map

```text
include/chargefw/  Public C++ API
src/core/          Toolkit-neutral molecular graph and conformers
src/features/      Derived topology, geometry, and spatial fragments
src/parameters/    Parameter models, JSON loading, and classification
src/methods/       Method interface, registry, applicability, and built-ins
src/calculation/   Selection, execution policy, full/cutoff dispatch, facade
src/charges/       Charge result types and validation
src/adapters/      Native and Gemmi-backed import/export adapters
apps/chargefw/     Command-line application
tests/cpp/         CTest executables and test support
cmake/             Dependencies, warnings, sanitizers, and install rules
old/               Archived ChargeFW2 compatibility reference
```

`include/chargefw/adapters/{native,gemmi}/all.h` are application convenience headers. Library code
should include only the individual adapter headers it uses.

## Architectural rules

- Keep `core` small and toolkit-neutral. File formats, toolkit objects, preparation policy, and cached
  derived data do not belong in `core::Molecule`.
- Put cached topology/geometry data and spatial fragments in `features`.
- Parameter matching is immutable external data:
  `ParameterSet + Molecule + TopologyFeatures -> ParameterClassification -> ParameterView`.
- Methods are stateless algorithms. Request-specific options, geometry, targets, and parameters enter
  through `methods::CalculationInput`.
- Keep scientific applicability, execution availability, deterministic selection, and calculation
  distinct, while allowing the application facade to compose them.
- Keep execution policy separate from method options. Never silently change method, parameter set,
  classification, execution mode/radius, charge correction, topology, protonation, or geometry.
- Preserve source atom order, molecule/conformer identity, and mappings at every result boundary.
- Missing scientific prerequisites are hard failures. Resource thresholds only guide automatic
  execution; explicit full execution may override them with a reported warning.
- Full calculations must remain faithful to their cited methods. Approximate execution must be
  capability-checked, explicit in provenance, and validated per method.
- Bindings, adapters, and applications translate representations into the native facade; they must not
  duplicate applicability, selection, or scientific policy.
- Avoid global mutable state. The library must remain suitable for concurrent and embedded use.

## ChargeFW2 compatibility

`old/` is read-only research material unless a task explicitly requests a compatibility fixture or
reference correction. Compare migrated behavior using identical graphs, formal charges, coordinates,
options, and parameter data. Record tolerances and intentional deviations in tests or `PROJECT.md`.
Do not merge prototype or archived behavior wholesale.

## Change discipline

1. Inspect the relevant public headers, implementation, tests, and `PROJECT.md` section first.
2. Make the smallest coherent change and keep unrelated formatting out of the diff.
3. Add focused tests for behavior changes, especially numerical, mapping, applicability, and policy
   cases.
4. Build the directly affected target, run its focused CTest test, then run `ctest --preset gcc-debug`.
   Use `clang-asan` for ownership/mapping changes and before completing reduced-execution work.
5. Stop for approval before heavyweight dependencies, broad rewrites, destructive commands, automatic
   chemistry policy, or public API changes beyond the requested scope.
6. Remove a TODO only when its full deliverable and validation are complete. Put durable implemented
   state in `PROJECT.md`.

## C++ conventions

- Use C++23. Public names remain under
  `chargefw::{core,features,parameters,methods,charges,calculation,adapters}`.
- Follow `.clang-format`: LLVM-derived, four spaces, 100-column limit, left-attached pointers and
  references.
- Prefer value ownership, `std::span` for read-only views, and explicit non-owning references/pointers
  when caller lifetime is guaranteed.
- Use exceptions for invalid inputs and calculation failures; preserve actionable target context.
