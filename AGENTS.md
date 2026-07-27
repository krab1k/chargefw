# ChargeFW Development Guide

## Mission

ChargeFW is a modern C++23 framework for empirical partial atomic-charge calculation. It is a
clean, library-first successor to the application-oriented ChargeFW2 engine used by Atomic Charge
Calculator III (ACC III). Preserve published method definitions and make preparation, parameter
selection, calculation, diagnostics, and integration explicit and reproducible.

The core must remain independent of file formats, RDKit, Gemmi, Python, GUI frameworks, web
services, and MCP. Those belong in optional adapters or applications.

## Repository map

```text
include/chargefw/  Public library API
src/core/          Molecule, atom, bond, conformer, periodic table
src/features/      Derived topology and conformer features
src/parameters/    Parameter models, JSON I/O, classification
src/methods/       Method interface, registry, applicability, built-ins
src/charges/       Charge result types and validation
apps/chargefw/     CLI application (currently a water demonstration)
tests/cpp/         CTest executables and test support
data/parameters/   Bundled empirical parameter sets
cmake/             Dependencies, tooling, warnings, sanitizers, install rules
old/               Archived ChargeFW2 reference implementation; do not modify casually
PROJECT.md         Current state, migration goals, development roadmap, and commands
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
- Keep applicability explicit. New requirements, limitations, and failures must be represented by
  `MethodRequirements`, prerequisite issues, and `ApplicabilityResult`, not hidden fallback logic.
- Preserve atom ordering and conformer identity in every calculation and output boundary.
- Never silently change protonation, geometry, parameter set, method, cutoff/cover policy, or
  charge correction. Report the selected policy and diagnostics.
- Keep method implementations faithful to their cited publications. Framework changes must not
  quietly change the reference/full calculation.

## ChargeFW2 compatibility

`old/` is the practical compatibility reference and contains 22 methods. The new registry
currently contains 19. Before claiming method parity, implement and validate:

- `sqe` (Split-charge Equilibration);
- `sqeq0` (SQE with initial formal charges);
- `sqeqp` (SQE with parameterized initial charges);
- their nine archived parameter sets.

For every migrated method, compare new and old implementations using identical molecule graphs,
formal charges, coordinates, options, and parameter files. Record intentional deviations and
numerical tolerances in tests or dedicated compatibility fixtures.

ChargeFW2 has application-coupled `full`, `cutoff`, and `cover` execution modes. Reproduce and
test compatible behavior before replacing it with a more general execution-policy abstraction.

## Change discipline

1. Read `PROJECT.md` and relevant public headers before designing a change.
2. Prefer the smallest coherent change; keep unrelated formatting out of diffs.
3. Add or update focused tests for behavior changes, especially numerical and applicability cases.
4. Validate the narrowest target first, then run the relevant CTest preset.
5. Do not add heavyweight dependencies, public APIs, broad rewrites, destructive commands, or
   automatic chemistry policies without explicit approval.
6. Do not modify `old/` except for an explicitly requested compatibility fixture or reference fix.

## C++ conventions

- C++23; public headers live in `include/chargefw/`.
- Follow `.clang-format`: LLVM-derived, four spaces, 100-column limit, left-attached pointers and
  references.
- Use value ownership for domain data, `std::span` for read-only views, and explicit non-owning
  references/pointers where lifetime is guaranteed by the caller.
- Use exceptions for invalid inputs and calculation failures; preserve actionable context in error
  messages.
- Keep public type names under `chargefw::{core,features,parameters,methods,charges}`.
- Avoid global mutable state. This library is intended for concurrent and embedded use.

## Build and test commands

Use CMake presets from the repository root:

```bash
# Debug build and full test suite
cmake --preset gcc-debug
cmake --build --preset gcc-debug
ctest --preset gcc-debug

# Sanitizer build and tests
cmake --preset clang-asan
cmake --build --preset clang-asan
ctest --preset clang-asan

# Static analysis build (requires clang-tidy)
cmake --preset clang-tidy
cmake --build --preset clang-tidy
```

See `PROJECT.md` for the full command reference, current known limitations, and development
priorities.
