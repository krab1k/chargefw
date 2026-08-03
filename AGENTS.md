# ChargeFW Development Guide

## Required reading and document ownership

Read this file before making changes. Then read the relevant sections of:

1. [PROJECT.md](PROJECT.md) for the current architecture, method coverage, compatibility goals,
   and product direction;
2. [TODO.md](TODO.md) for actionable, categorized deliverables and completion criteria;
3. [README.md](README.md) for the supported build, test, installation, and formatting commands.

Keep these documents distinct: this file defines implementation rules; `PROJECT.md` defines the
technical/project state; `TODO.md` tracks work; and `README.md` is the concise entry point. Update
the owning document when a change makes its content stale rather than duplicating it elsewhere.

When implementation changes architecture, supported methods, compatibility status, build/run/test
commands, installation behavior, public workflow, or planned deliverables, update the relevant
owned document in the same change. Do not leave documentation updates as an implicit follow-up.

## Repository map

```text
include/chargefw/  Public library API
                   (`adapters/{native,gemmi}/all.h` are per-backend convenience umbrella
                   headers for applications; library code should continue to include only the
                   individual headers it uses)
src/core/          Molecule, atom, bond, conformer, periodic table
src/features/      Derived topology and conformer features
src/calculation/   High-level applicability, selection, and execution facade
src/parameters/    Parameter models, JSON I/O, classification
src/methods/       Method interface, registry, applicability, built-ins
src/charges/       Charge result types and validation
apps/chargefw/     CLI application (currently a water demonstration)
tests/cpp/         CTest executables and test support
cmake/             Dependencies, tooling, warnings, sanitizers, install rules
old/               Archived ChargeFW2 reference implementation; do not modify casually
```

## Architectural rules

- Keep `core` small and toolkit-neutral. `core::Molecule` is a graph with optional conformers,
  not an SDF, RDKit, PDB, or mmCIF object.
- Put cached derived data in `features` (`TopologyFeatures`, `ConformerFeatures`), not in the
  core molecule model.
- Treat parameter matching as immutable external data:
  `ParameterSet + Molecule + TopologyFeatures -> ParameterClassification -> ParameterView`.
  Do not mutate molecules while evaluating parameter sets.
- Methods are stateless algorithms. Pass options, geometry, and parameter data through
  `methods::CalculationInput`; do not store mutable request-specific state in a `Method`.
- Keep applicability explicit. Model requirements, limitations, and failures with
  `MethodRequirements`, prerequisite issues, and `ApplicabilityResult`; do not hide fallback
  behavior.
- Preserve atom ordering and conformer identity at every calculation and output boundary.
- Never silently change protonation, geometry, parameter set, method, cutoff/cover policy, or
  charge correction. Report the selected policy and diagnostics.
- Keep methods faithful to their cited publications. Framework changes must not quietly alter the
  reference/full calculation.
- Keep the core independent of file formats, RDKit, Gemmi, Python, GUI frameworks, web services,
  and MCP. Put those concerns in optional adapters or applications.
- Make bindings and adapters translate representations into the native calculation facade; do not
  duplicate applicability, selection, parameter, or scientific calculation policy in integrations.

## ChargeFW2 compatibility

`old/` is the compatibility reference. Do not modify it except for an explicitly requested
compatibility fixture or reference fix. For every migrated method, compare the new and archived
implementations with identical graphs, formal charges, coordinates, options, and parameter files.
Record intentional deviations and numerical tolerances in tests or dedicated compatibility
fixtures. See [PROJECT.md](PROJECT.md#compatibility-gap-with-chargefw2) and
[TODO.md](TODO.md#method-parity-and-scientific-validation) for the current parity scope.

## Change discipline

1. Read relevant public headers and the applicable `PROJECT.md` section before designing a change.
2. Prefer the smallest coherent change; keep unrelated formatting out of diffs.
3. Add or update focused tests for behavior changes, especially numerical and applicability cases.
4. Run the narrowest relevant target first, then the applicable CTest preset from `README.md`.
5. Do not add heavyweight dependencies, public APIs, broad rewrites, destructive commands, or
   automatic chemistry policies without explicit approval.
6. Keep implementation changes and TODO completion synchronized: check an item only when its
   stated deliverable and validation are complete.

## C++ conventions

- Use C++23. Public headers live in `include/chargefw/` and public type names remain under
  `chargefw::{core,features,parameters,methods,charges,calculation}`.
- Follow `.clang-format`: LLVM-derived, four spaces, 100-column limit, left-attached pointers and
  references.
- Prefer value ownership for domain data, `std::span` for read-only views, and explicit
  non-owning references/pointers when the caller guarantees lifetime.
- Use exceptions for invalid inputs and calculation failures; preserve actionable error context.
- Avoid global mutable state. The library must be suitable for concurrent and embedded use.
