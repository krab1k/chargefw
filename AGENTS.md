# ChargeFW Development Guide

## Required reading and document ownership

Read this file before changing the repository. Then read only the relevant sections of:

- [docs/PROJECT.md](docs/PROJECT.md): implemented architecture, philosophy, capabilities, and limits;
- [docs/CLI.md](docs/CLI.md): implemented command-line behavior;
- [docs/NATIVE.md](docs/NATIVE.md): implemented public C++ API and installation;
- [docs/PYTHON.md](docs/PYTHON.md): implemented Python API and package behavior;
- [TODO.md](TODO.md): unfinished product work; and
- [README.md](README.md): concise project entry point.

Keep ownership clear. The files under `docs/` are user documentation and must describe implemented
behavior, not development plans or milestone history. Unfinished work belongs in `TODO.md`; remove it
when complete rather than retaining checked history. `README.md` should remain a concise user entry
point.

## Repository map

```text
include/chargefw/   Public C++ API
src/                Native library implementation, organized by public domain
apps/chargefw/      Command-line application
python/chargefw/    Public Python package and private extension stubs
python/src/         Nanobind extension implementation
data/parameters/    Bundled parameter sets
tests/cpp/          Native tests and downstream CMake consumer
tests/python/       Python package and adapter tests
tests/fixtures/     Molecular and parameter test data
cmake/              Dependencies, diagnostics, and installation rules
docs/               Implemented user documentation
Dockerfile           User-facing CLI container image
docker/              Developer compatibility containers
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

## Change discipline

1. Inspect the relevant public headers, implementation, tests, and user documentation section first.
2. Make the smallest coherent change and run clang-format each time on modified files.
3. Add focused tests for behavior changes, especially numerical, mapping, applicability, and policy
   cases.
4. Follow the validation cadence below. Build the directly affected target and run its focused test first;
   do not make the full compiler/sanitizer matrix part of every edit cycle.
5. Stop for approval before heavyweight dependencies, broad rewrites, destructive commands, automatic
   chemistry policy, or public API changes beyond the requested scope.
6. Remove a TODO only when its full deliverable and validation are complete. Update the owning user
   document with durable implemented behavior.

### Validation cadence

- Use `gcc-debug` for the rapid edit/build/focused-test loop. After a coherent change, run the full
  `ctest --preset gcc-debug` suite.
- Use `clang-debug` for cross-compiler validation when changing templates, conversions, overloads,
  headers, or compiler-sensitive C++ behavior; run its focused test before substantial work is complete.
- Use `gcc-release` and `clang-release` for optimization- or `NDEBUG`-sensitive behavior, numerical
  methods, Eigen code, and before completing substantial changes. Run the affected test under both
  release configurations when practical.
- Use `clang-asan` for ownership, lifetime, bounds, mapping, parser, container, view, and pointer
  changes. Use `clang-ubsan` for arithmetic, conversions, shifts, alignment, indexing, and other
  undefined-behavior risks. Run the focused sanitizer test first and the full relevant sanitizer suite
  before completing substantial risk-sensitive work. Keep sanitizer builds serial or deliberately
  low-parallel because template-heavy translation units can consume substantial memory.
- Run `clang-tidy` after meaningful implementation or public-interface changes and before a milestone;
  it is static analysis, not a replacement for compiler or runtime tests.
- Before a substantial merge or milestone, run the full debug, release, sanitizer, and clang-tidy
  matrix. CI may distribute those configurations across independent jobs.

### Container validation

The root `Dockerfile` builds the user-facing Ubuntu 26.04 CLI image. It performs an optimized native
release build without tests, Python bindings, or host-specific `-march=native` instructions.

`docker/Dockerfile.ubuntu-26.04`, `docker/Dockerfile.debian-13`, and `docker/Dockerfile.fedora-44` are
developer compatibility containers. Each builds and tests the `gcc-release` configuration; run them one
at a time because template-heavy builds consume substantial memory. Fedora is pinned to release 44; update
that release intentionally before it reaches end of life.

## C++ conventions

- Use C++23. Public names remain under
  `chargefw::{core,features,parameters,methods,charges,calculation,adapters}`.
- Follow `.clang-format`: LLVM-derived, four spaces, 100-column limit, left-attached pointers and
  references.
- Prefer value ownership, `std::span` for read-only views, and explicit non-owning references/pointers
  when caller lifetime is guaranteed.
- Use exceptions for invalid inputs and calculation failures; preserve actionable target context.
