# ChargeFW

ChargeFW is a C++23, library-first framework for empirical partial atomic-charge calculation.
It is a modern successor to ChargeFW2, the engine used by Atomic Charge Calculator III (ACC III).

The project currently provides a toolkit-neutral molecular core, built-in empirical methods,
parameter-set loading and classification, applicability checks, structured charge results, and
native MOL/SDF/MOL2 and JSON input adapters. The `chargefw` executable is a small file demonstration; it does not
yet provide a full user-facing file/SMILES CLI.

## Documentation map

Read the document appropriate to the task instead of duplicating information between them:

| Document | Audience | Source of truth for |
| --- | --- | --- |
| [README.md](README.md) | Users and contributors | Prerequisites, first build, demo, and tests. |
| [PROJECT.md](PROJECT.md) | Maintainers and designers | Current architecture, method coverage, compatibility context, and product roadmap. |
| [AGENTS.md](AGENTS.md) | AI agents and implementers | Repository boundaries, architectural constraints, coding rules, and change discipline. |
| [TODO.md](TODO.md) | Maintainers and planners | Categorized, actionable deliverables toward a production-quality charge-calculation tool. |

For implementation work, read [AGENTS.md](AGENTS.md) first, then the relevant section of
[PROJECT.md](PROJECT.md), and use [TODO.md](TODO.md) to identify scoped deliverables.

## Requirements

- CMake 3.27 or newer
- Ninja
- GCC or Clang with C++23 support
- Internet access on the first configure when Eigen 5.0.1 and nlohmann/json 3.12.0 are not
  already available to CMake

CMake first looks for Eigen3 and nlohmann/json on the system, then obtains them with
`FetchContent` if needed.

## Build, run, and test

Run commands from the repository root:

```bash
# Configure and build the debug configuration.
cmake --preset gcc-debug
cmake --build --preset gcc-debug

# Run the complete debug test suite.
ctest --preset gcc-debug

# Run the molecular-file demonstration from the build tree.
CHARGEFW_PARAMETER_DIR="$PWD/data/parameters" build/gcc-debug/apps/chargefw/chargefw tests/water.sdf
```

Run an individual test after building:

```bash
ctest --test-dir build/gcc-debug -R test_qeq --output-on-failure
```

The demo accepts `.sdf`, `.mol`, `.mol2`, and ChargeFW `.json` input files, selects the reader from
the file extension, reports and skips malformed records, loads bundled parameter sets, and
autodetects the highest-priority applicable method and parameter set. It writes a versioned JSON
result document to standard output, or to a file passed with `--output`.

The JSON result uses `"schema_version": "1.0"` and retains one `results` entry for every successfully
imported molecule. Each entry records input identity, explicit atom and conformer mappings, selected
method and optional parameter-set IDs, atom-order-preserving charge assignments, and structured
diagnostics. Native input adapters currently preserve atom and conformer order, so their mappings
are reported as `{ "kind": "identity" }`.

JSON partial-charge values are rounded to at most four decimal places; calculations retain their
native precision internally.

Native adapter headers are explicitly directional: `json_input`, `mol_input`, `sdf_input`, and
`mol2_input` read molecule records, while `json_output::JsonWriter` writes calculation-result
documents.

ChargeFW JSON input uses `"schema_version": "1.0"` and a `molecules` array. Each molecule has
required `atoms` with `atomic_number` and `formal_charge`, optional indexed `bonds`, and optional
conformers with coordinate triplets. Array order is preserved; bonds are never inferred from geometry.

## Other configurations

```bash
# Optimized build; tests are disabled by this preset.
cmake --preset gcc-release
cmake --build --preset gcc-release

# AddressSanitizer + UndefinedBehaviorSanitizer build and tests.
cmake --preset clang-asan
cmake --build --preset clang-asan
ctest --preset clang-asan

# Local installation to _install, then test without an environment override.
cmake --preset local-install
cmake --build build/local-install
cmake --install build/local-install
env -u CHARGEFW_PARAMETER_DIR _install/bin/chargefw tests/water.sdf
```

The `_install` directory is a local installation staging area for testing installation behavior; it
is not intended as a system-wide install location. The installed executable locates its parameter
directory automatically. The build-tree executable requires `CHARGEFW_PARAMETER_DIR` as shown
above. For the complete product plan, see [PROJECT.md](PROJECT.md).

Configure options can be passed to a preset, for example:

```bash
cmake --preset gcc-debug -DCHARGEFW_BUILD_TESTS=OFF
cmake --preset gcc-debug -DCHARGEFW_BUILD_CLI=OFF
cmake --preset gcc-debug -DCHARGEFW_ENABLE_CCACHE=OFF
cmake --preset clang-debug -DCHARGEFW_ENABLE_SANITIZERS=ON
cmake --preset clang-debug -DCHARGEFW_ENABLE_CLANG_TIDY=ON
```

## Formatting

The optional [pre-commit](https://pre-commit.com/) hook formats staged C++ files using the
repository `.clang-format` configuration:

```bash
pre-commit install
```

Run `pre-commit run --all-files` to format all applicable files. Commit
`.pre-commit-config.yaml`, not the generated `.git/hooks/pre-commit` file.
