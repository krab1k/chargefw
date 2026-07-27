# ChargeFW

ChargeFW is a C++23 framework for empirical partial atomic-charge calculation. It is currently a
library-first project with a demonstration executable; it is not yet a general-purpose command-line
tool for molecular files or SMILES.

The demo constructs two water conformers in code, loads the bundled parameter sets, finds
applicable methods, and prints calculated atomic charges.

For architecture, method coverage, and development plans, see [PROJECT.md](PROJECT.md).

## Requirements

- CMake 3.27 or newer
- Ninja
- GCC or Clang with C++23 support
- Internet access on the first configure when Eigen 5.0.1 and nlohmann/json 3.12.0 are not
  already available to CMake

CMake first looks for Eigen3 and nlohmann/json on the system, then obtains them with
`FetchContent` if needed.

## Configure and build

Run these commands from the repository root:

```bash
cmake --preset gcc-debug
cmake --build --preset gcc-debug
```

This creates a debug build in `build/gcc-debug`. Tests and the demonstration executable are
enabled by default.

To build a release configuration:

```bash
cmake --preset gcc-release
cmake --build --preset gcc-release
```

The release preset writes to `build/gcc-release` and disables tests.

## Run the demonstration

The build-tree executable needs the location of the bundled parameter data:

```bash
CHARGEFW_PARAMETER_DIR="$PWD/data/parameters" build/gcc-debug/apps/chargefw/chargefw
```

The program reports the loaded methods and parameter sets, then prints charge assignments for the
applicable methods. Its input is fixed water data compiled into the executable; command-line input
files and SMILES are not supported yet.

## Install and run

The local preset installs to `_install` in the repository:

```bash
cmake --preset clion-local
cmake --build build/clion-local
cmake --install build/clion-local
_install/bin/chargefw
```

The installed executable finds the installed parameter directory automatically.

To choose a different install location, configure with a custom `CMAKE_INSTALL_PREFIX`:

```bash
cmake --preset gcc-release -DCMAKE_INSTALL_PREFIX=/path/to/install
cmake --build --preset gcc-release
cmake --install build/gcc-release
```

## Test

```bash
ctest --preset gcc-debug
```

Run one test by name after building:

```bash
ctest --test-dir build/gcc-debug -R test_qeq --output-on-failure
```

## Automatic C++ formatting before commit

The repository includes a [pre-commit](https://pre-commit.com/) hook that formats staged C++ files
with `clang-format` and the repository's `.clang-format` configuration. Install `pre-commit` and
`clang-format`, then enable the hook once per clone:

```bash
pre-commit install
```

Future `git commit` commands automatically format the affected C++ files. Run it manually across
the repository with:

```bash
pre-commit run --all-files
```

Commit `.pre-commit-config.yaml` so all contributors use the same formatting hook; the generated
`.git/hooks/pre-commit` file remains local.

## Useful configuration options

Pass options when configuring a preset:

```bash
cmake --preset gcc-debug -DCHARGEFW_BUILD_TESTS=OFF
cmake --preset gcc-debug -DCHARGEFW_BUILD_CLI=OFF
cmake --preset gcc-debug -DCHARGEFW_ENABLE_CCACHE=OFF
cmake --preset clang-debug -DCHARGEFW_ENABLE_SANITIZERS=ON
cmake --preset clang-debug -DCHARGEFW_ENABLE_CLANG_TIDY=ON
```

Available project options:

| Option | Default | Purpose |
| --- | --- | --- |
| `CHARGEFW_BUILD_TESTS` | `ON` | Build the CTest suite. |
| `CHARGEFW_BUILD_CLI` | `ON` | Build the `chargefw` demonstration executable. |
| `CHARGEFW_ENABLE_CCACHE` | `ON` | Use `ccache` when it is available. |
| `CHARGEFW_ENABLE_WARNINGS` | `ON` | Enable project compiler warnings. |
| `CHARGEFW_ENABLE_SANITIZERS` | `OFF` | Enable AddressSanitizer and UndefinedBehaviorSanitizer. |
| `CHARGEFW_ENABLE_CLANG_TIDY` | `OFF` | Run clang-tidy during compilation. |

## Additional build presets

```bash
# Clang debug build
cmake --preset clang-debug
cmake --build --preset clang-debug

# Clang build with AddressSanitizer and UndefinedBehaviorSanitizer
cmake --preset clang-asan
cmake --build --preset clang-asan
ctest --preset clang-asan

# Run clang-tidy during compilation (clang-tidy must be installed)
cmake --preset clang-tidy
cmake --build --preset clang-tidy
```
