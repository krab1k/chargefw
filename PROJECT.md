# ChargeFW Project Guide

## Purpose

ChargeFW is a C++23 framework for empirical partial atomic-charge calculation. It provides a
toolkit-neutral core model, built-in empirical methods, parameter-set loading and classification,
method applicability checks, and structured charge results.

It is being developed as a modern, library-first successor to **ChargeFW2**, the current
computational engine of [Atomic Charge Calculator III (ACC III)](https://acc.biodata.ceitec.cz).
ACC III is described in Raček *et al.*, *Atomic Charge Calculator III: a modern platform for
calculating partial atomic charges*, **Nucleic Acids Research** (2026),
DOI [10.1093/nar/gkag379](https://doi.org/10.1093/nar/gkag379).

ChargeFW is not yet the ACC III backend. `old/` is an archived ChargeFW2 copy used for behavior,
parameter, and compatibility research.

## Current architecture

```text
External molecular adapters / applications                 Not yet implemented
    (SDF, Mol2, PDB, mmCIF, SMILES, RDKit, Gemmi)
                              |
                              v
core::MoleculeCollection
    |- core::Molecule: atoms, bonds, conformers
    |- core::Atom: atomic number, formal charge, source name
    |- core::Bond: atom indices and order
    `- core::Conformer: coordinates sharing the molecule topology
                              |
                              v
features::PreparedMoleculeCollection
    |- PreparedMolecule: molecule + cached TopologyFeatures
    `- ConformerFeatures: on-demand geometry view per conformer
                              |
                 +------------+------------+
                 |                         |
                 v                         v
     parameters::ParameterSet       methods::MethodRegistry
                 |                         |
                 v                         v
      ParameterClassification     MethodRequirements/options
                 |                         |
                 +------------+------------+
                              v
                methods::find_applicable_methods()
                              |
                              v
                      methods::ApplicableMethod
                              |
                              v
                   methods::calculate_charges()
                              |
                              v
                      charges::ChargeSet
```

### Important public types

| Area | Types | Role |
|---|---|---|
| Core | `Atom`, `Bond`, `Conformer`, `Molecule`, `MoleculeCollection` | Input molecular graph and coordinates. |
| Features | `TopologyFeatures`, `ConformerFeatures`, `PreparedMolecule` | Cached/derived topology and geometry. |
| Parameters | `ParameterSet`, `ParameterClassification`, `ParameterView` | Parameter storage, matching, and method-facing lookup. |
| Methods | `Method`, `MethodRegistry`, `MethodRequirements`, `MethodOptions`, `ApplicableMethod` | Algorithm interface, capabilities, selection, and execution. |
| Charges | `AtomicCharges`, `ChargeAssignment`, `ChargeSet`, `ChargeCollection` | Atom-indexed calculated results and provenance. |

## Typical library workflow

```cpp
core::MoleculeCollection molecules{/* validated molecules */};
features::PreparedMoleculeCollection prepared{molecules};

const auto parameter_sets = parameters::load_default_parameter_sets();
const auto& registry = methods::method_registry();

std::vector<const methods::Method*> candidates;
for (const auto& method : registry.methods()) {
    candidates.push_back(method.get());
}

const auto applicability =
    methods::find_applicable_methods(prepared, candidates, parameter_sets);

// Application policy selects an explicitly reported candidate.
const auto& selected = applicability.applicable.front();
const charges::ChargeSet result = methods::calculate_charges(selected, prepared);
```

The current `chargefw` executable is a **water demonstration**. It builds two water conformers in
code, loads bundled parameter sets, identifies applicable methods, chooses the lowest numeric
parameter priority per method, and prints charges. It is not yet a user-facing file/SMILES CLI.

## Implemented methods and parameters

The current registry contains 19 methods:

```text
abeem, charge2, delre, denr, dummy, eem, eqeq, eqeqc, formal, gdac,
kcm, mgc, mpeoe, peoe, qeq, sfkeem, smpqeq, tsef, veem
```

Bundled JSON parameter sets cover these parameterized methods and variants. `data/parameters/`
is installed under `share/chargefw/parameters`.

### Compatibility gap with ChargeFW2

The archived ChargeFW2 registry contains the same 19 methods plus:

```text
sqe, sqeq0, sqeqp
```

It also contains nine corresponding SQE-family parameter files. Method and parameter parity is a
release-blocking compatibility objective because ACC III highlights SQE+qp.

## ChargeFW2 research summary

ChargeFW2 is a functional, application-oriented engine with CLI/Python bindings, custom
SDF/Mol2 reading, Gemmi PDB/mmCIF reading, output writers, OpenMP, and a nanoflann KD-tree. It
has useful production behavior, but combines concerns that are intentionally separated here:

| ChargeFW2 pattern | ChargeFW replacement |
|---|---|
| Atom stores graph data, coordinates, PDB residue data, Mol2 types, and classification state. | Core graph stays minimal; adapters and derived feature layers own source-specific metadata/caches. |
| Molecule owns all-pairs topology matrices and spatial KD-tree. | `TopologyFeatures` and `ConformerFeatures` own derived data. |
| Parameter classification mutates atoms and bonds. | Immutable `ParameterClassification` and `ParameterView`. |
| Method owns mutable parameters pointer and option values. | Stateless `Method::calculate(CalculationInput)`. |
| Suitability is mainly filtering/console behavior. | Structured applicability and prerequisite diagnostics. |
| PDB/mmCIF reader uses the first model only. | Core can represent multiple conformers; adapters should define model/altloc policy. |
| `full`, `cutoff`, and `cover` are coupled inside `EEMethod`. | Future execution policy should be explicit, reusable, and provenance-rich. |

The new implementation must preserve validated scientific behavior before making intentional
improvements. Compatibility fixtures should use identical topology, coordinates, formal charge,
options, and parameter data.

## Development priorities

### Near-term: credible core release

1. Implement SQE, SQE+q0, and SQE+qp with archived parameter-set parity.
2. Build a ChargeFW2-to-ChargeFW numerical compatibility suite for every method and parameter set.
3. Document numerical tolerances and intentional behavioral deviations.
4. Add a high-level calculation request/result facade and selection policy outside the demo CLI.
5. Add CMake package export (`find_package(chargefw CONFIG REQUIRED)`).
6. Replace the demo with a real CLI and JSON result/provenance output.

### Integration layer, driven by actual users

1. Optional RDKit adapter: SMILES, SDF/MOL V2000/V3000, Mol2, hydrogen/conformer preparation.
2. SDF batch streaming with bounded memory, record-level errors, and deterministic output order.
3. Optional Gemmi adapter: PDB/mmCIF with explicit component/model/altloc/bond policies.
4. Python bindings once the native request/result boundary is stable.
5. ACC III shadow comparison, then method-by-method adoption only after parity is established.

### Future work; do not pre-build

- cutoff/cover redesign after parity, benchmarks, and an explicit execution-policy contract;
- WebAssembly JSON API only if a browser consumer needs local calculation;
- MCP server only as an adapter over stable JSON/library APIs;
- ML/GNN model adapters as optional estimators, not a dependency of `chargefw_core`;
- custom parsers only where optional RDKit/Gemmi adapters are unsuitable.

## Build requirements

- CMake 3.27 or newer;
- Ninja;
- GCC or Clang with C++23 support;
- internet access on first configure if Eigen 5.0.1 or nlohmann/json 3.12.0 are not installed;
- optional: `ccache`, `clang-tidy`.

CMake finds `Eigen3` and `nlohmann_json` first, then downloads them with `FetchContent` when
absent. Tests and CLI are enabled by default.

## CMake presets and commands

All commands are run from the repository root.

### Debug build and tests

```bash
cmake --preset gcc-debug
cmake --build --preset gcc-debug
ctest --preset gcc-debug
```

The configured build directory is `build/gcc-debug`. To run an individual test after building:

```bash
ctest --test-dir build/gcc-debug -R test_qeq --output-on-failure
```

### Release build

```bash
cmake --preset gcc-release
cmake --build --preset gcc-release
```

This preset disables tests and writes to `build/gcc-release`.

### Clang, sanitizers, and clang-tidy

```bash
# Clang debug build
cmake --preset clang-debug
cmake --build --preset clang-debug

# AddressSanitizer + UndefinedBehaviorSanitizer and test suite
cmake --preset clang-asan
cmake --build --preset clang-asan
ctest --preset clang-asan

# clang-tidy during compilation; clang-tidy must be installed
cmake --preset clang-tidy
cmake --build --preset clang-tidy
```

### Local install and demo executable

`CMakeUserPresets.json` provides a local install preset:

```bash
cmake --preset clion-local
cmake --build build/clion-local
cmake --install build/clion-local

# Run from the build tree
CHARGEFW_PARAMETER_DIR="$PWD/data/parameters" build/clion-local/apps/chargefw/chargefw

# Or after installation
_install/bin/chargefw
```

The build-tree executable needs `CHARGEFW_PARAMETER_DIR` because bundled data is installed rather
than copied beside it. The installed executable locates its installed parameter directory
automatically. The existing `_install/` directory is local development output and is ignored by
Git.

### Configure options

Useful cache options:

```bash
cmake --preset gcc-debug -DCHARGEFW_BUILD_TESTS=OFF
cmake --preset gcc-debug -DCHARGEFW_BUILD_CLI=OFF
cmake --preset gcc-debug -DCHARGEFW_ENABLE_CCACHE=OFF
cmake --preset clang-debug -DCHARGEFW_ENABLE_SANITIZERS=ON
cmake --preset clang-debug -DCHARGEFW_ENABLE_CLANG_TIDY=ON
```

## Validation expectations

For any scientific-method change:

1. Add focused unit/regression tests.
2. Run the narrowest related test executable or CTest regex.
3. Run `ctest --preset gcc-debug` before considering the change complete.
4. Run `ctest --preset clang-asan` for memory-sensitive, feature-cache, parser/adapter, or
   numerical changes where practical.
5. For migrated ChargeFW2 behavior, add a parity fixture rather than relying on manual output
   comparison.

## Scientific and integration principles

- Preserve source atom indexing and map every output charge back to its source atom.
- Record method, parameter set, method options, conformer, approximation policy, and warnings in
  every serializable result.
- Keep exact/reference calculations distinct from cutoff/cover approximations.
- Make automatic selection explainable and overridable.
- Treat empirical charge results as method-dependent estimates; expose coverage and numerical
  limitations rather than implying universal validity.
