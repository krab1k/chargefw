# Native C++ library

ChargeFW exposes a C++23 shared library through the `chargefw::core` CMake target. The public API owns
molecules, parameter sets, assessments, plans, results, and charge arrays; ordinary application code does
not need to manage prepared-feature or classification lifetimes.

## Build and link

ChargeFW requires CMake 3.28 or newer and a GCC or Clang toolchain with C++23 support.

```bash
cmake --preset gcc-release -DCMAKE_INSTALL_PREFIX="$PWD/_install"
cmake --build --preset gcc-release
cmake --install build/gcc-release --strip
```

Use the installed CMake package from another project:

```cmake
find_package(chargefw CONFIG REQUIRED)
target_link_libraries(my_application PRIVATE chargefw::core)
```

When ChargeFW builds its bundled dependencies, the installation includes the public nlohmann/json and
Gemmi development packages and the private oneTBB runtime. With
`CHARGEFW_USE_SYSTEM_DEPENDENCIES=ON`, those dependency packages must remain discoverable by downstream
CMake.

The installed prefix can be relocated as a unit. `load_default_parameter_sets()` resolves bundled JSON
relative to the loaded ChargeFW library.

## Basic calculation

```cpp
#include <chargefw/calculation/calculation.h>
#include <chargefw/core/atom.h>
#include <chargefw/core/bond.h>
#include <chargefw/core/conformer.h>
#include <chargefw/core/molecule.h>
#include <chargefw/core/molecule_collection.h>
#include <chargefw/parameters/io/parameter_set_io.h>

#include <iostream>
#include <utility>
#include <vector>

int main() {
    using namespace chargefw;

    auto atoms = std::vector<core::Atom>{core::Atom{8}, core::Atom{1}, core::Atom{1}};
    auto bonds = std::vector<core::Bond>{core::Bond{0, 1}, core::Bond{0, 2}};
    auto conformers = std::vector<core::Conformer>{core::Conformer{{
        {0.0, 0.0, 0.0},
        {0.96, 0.0, 0.0},
        {-0.24, 0.93, 0.0},
    }}};
    auto molecule = core::Molecule{
        std::move(atoms), std::move(bonds), std::move(conformers), "water"};

    auto request = calculation::AssessmentRequest{
        .molecules = core::MoleculeCollection{
            std::vector<core::Molecule>{std::move(molecule)}},
        .parameter_sets = parameters::load_default_parameter_sets(),
        .method_id = "eem",
        .execution_selection = calculation::ExecutionSelection{
            calculation::ExecutionSelectionKind::full},
    };

    auto assessment = calculation::assess(std::move(request));
    auto result = calculation::calculate(assessment);
    if (!result.calculated()) {
        std::cerr << result.failure_message.value_or("calculation failed") << '\n';
        return 1;
    }

    for (double value : result.charges->assignment(0).charges.values()) {
        std::cout << value << '\n';
    }
}
```

Passing an lvalue `AssessmentRequest` copies its molecule and parameter inputs. Passing an rvalue
transfers them into the assessment. An `AssessmentResult` owns the prepared molecules, classifications,
and plans needed for execution.

## Molecule model

`chargefw::core::Molecule` owns:

- a source-ordered vector of `Atom` values with atomic number, formal charge, and optional source name;
- indexed single, double, or triple `Bond` values;
- zero or more `Conformer` values, each containing one position per atom; and
- an optional molecule name.

`MoleculeCollection` owns source-ordered molecules and an optional collection name. Constructors validate
atomic numbers, bond endpoints and duplicates, conformer cardinality, and finite coordinates. Methods that
require geometry are inapplicable to molecules without conformers.

## Assessment

`calculation::AssessmentRequest` contains:

- the owned `MoleculeCollection` and parameter sets;
- optional method and parameter-set IDs;
- method-scoped `methods::MethodOptions` overrides;
- strict or permissive parameter classification options;
- an `ExecutionSelection`; and
- a `ResourcePolicy` containing automatic thresholds and the default execution thread limit.

`calculation::assess()` returns an `AssessmentResult` with priority-ordered `plans()`, structured
`rejections()`, `default_plan()`, and applicability timing. A plan exposes its applicable candidate,
concrete `ExecutionPolicy`, and warnings. Plans are tied to the assessment that created them.

Use `methods::method_registry()` to inspect built-in method metadata and option schemas. Use
`parameters::load_parameter_set_json*()` to load custom native parameter data, or
`load_default_parameter_sets()` for the installed catalog. Parameter-set IDs in one assessment request
must be unique.

## Execution policy

`ExecutionSelectionKind` accepts automatic, full, cutoff, or cover selection. Effective plans always use
the concrete `ExecutionMode` values full, cutoff, or cover.

- Full execution requires no radius and the `none` charge-correction policy.
- Explicit cutoff and cover require a finite radius of at least 8 Å.
- Automatic reduced execution uses 12 Å.
- The default full-to-cutoff and cutoff-to-cover thresholds are 20,000 and 80,000 atoms.
- `std::nullopt` disables the corresponding resource threshold.
- A thread count of zero delegates scheduling to oneTBB.

Automatic policy does not turn missing scientific requirements into warnings. Explicit unsupported
execution produces no runnable plan, and explicit full execution above a threshold remains runnable with
a resource warning.

## Results, errors, and mapping

`calculation::ExecutionResult` reports an `ExecutionStatus`, optional `charges::ChargeSet`, rejections,
effective calculation provenance, optional failure text, and applicability/computation timings.

The status is one of `success`, `invalid_input_or_request`, `no_executable_plan`, `numerical_failure`, or
`cancelled`. Invalid native API inputs generally throw `std::invalid_argument`; the owned facade converts
calculation failures and cooperative cancellation into result statuses. The lower-level
`CalculationRequest` overload propagates calculation and cancellation exceptions.

A `ChargeSet` owns the method ID, optional parameter-set ID, and source-indexed assignments.
Geometry-dependent methods return one assignment per molecule conformer. Geometry-independent methods
return one assignment per molecule. Every `ChargeAssignment` carries a molecule index and an optional
conformer index; its charge vector follows source atom order.

## Reuse, progress, and concurrency

`calculate(assessment, plan)` executes a concrete plan without repeating preparation, classification, or
applicability checks. Assessments and plans can be reused, and independent calculations can run
concurrently.

Derive from `calculation::CalculationObserver` to receive computation, target, and reduced-fragment
progress or to request cooperative cancellation. The observer is non-owning and callbacks may run on
oneTBB workers, so implementations must be thread-safe and must not throw. No partial charge set is
returned on cancellation.

## Molecular adapters

Individual public headers are provided under `chargefw/adapters/native` and
`chargefw/adapters/gemmi`. The `all.h` headers are convenience includes for applications.

- `MolReader`, `SdfReader`, and `Mol2Reader` import native molecular formats.
- `JsonReader` imports ChargeFW molecule JSON 1.0.
- `PdbReader` and `MmcifReader` import Gemmi structures with explicit record, bond, and conformer policy.
- Native writers emit ChargeFW JSON, SDF, and MOL2.
- The Gemmi writer emits preservation-oriented or generated mmCIF charge data.

Readers return `ImportedMoleculeRecord`, which keeps the molecule together with source identity and import
diagnostics. SDF and MOL2 readers consume one record at a time; the CLI currently materializes the full
collection before assessment.
