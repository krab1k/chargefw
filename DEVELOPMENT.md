# ChargeFW development

This is the practical guide to building and testing ChargeFW. Repository policy and architecture rules
are in [AGENTS.md](AGENTS.md); public behavior is documented under [`docs/`](docs/).

## Requirements

Native development requires CMake 3.28 or newer, Ninja, and a C++23-capable GCC or Clang toolchain.
CMake downloads pinned dependencies when compatible system packages are unavailable.

The checked-in presets also build the Python bindings. Their Python interpreter needs Python 3.10 or
newer, development headers, NumPy, Nanobind, mypy, and Gemmi 0.7.4. RDKit is optional. The presets use
`/usr/bin/python3`; select another interpreter explicitly when needed:

```bash
cmake --preset gcc-debug \
    -DCHARGEFW_PYTHON_EXECUTABLE=/path/to/python
```

## Build and test

Use `gcc-debug` for normal development:

```bash
cmake --preset gcc-debug
cmake --build build/gcc-debug
ctest --test-dir build/gcc-debug --output-on-failure -E '^cpptest$'
```

Gemmi registers its internal `cpptest` unconditionally even though that executable is excluded from
normal builds. Exclude it when invoking CTest directly.

For a quick edit cycle, build the affected target and run its test by name:

```bash
cmake --build build/gcc-debug --target test_molecule
ctest --test-dir build/gcc-debug --output-on-failure -E '^cpptest$' -R '^test_molecule$'
```

Available presets are:

| Preset | Use |
| --- | --- |
| `gcc-debug` | Primary development and test loop |
| `clang-debug` | Cross-compiler validation |
| `gcc-release`, `clang-release` | Optimized and numerical validation |
| `clang-asan` | Memory and lifetime errors |
| `clang-ubsan` | Undefined behavior |
| `clang-tidy` | Static analysis during compilation |

All checked-in presets build the native library, CLI, tests, and Python bindings. Release presets also
enable IPO and host-native optimization; debug, sanitizer, and static-analysis presets disable
host-native optimization for reproducibility. Ignored `CMakeUserPresets.json` presets may inherit the
stable `gcc-debug` and `gcc-release` names to customize local build and installation directories.

Configure another preset in the same way, then build and test its directory. `clang-tidy` runs during
its build and does not need a test run:

```bash
cmake --preset clang-tidy
cmake --build build/clang-tidy
```

Build sanitizer presets with low parallelism because they can use substantial memory:

```bash
cmake --preset clang-asan
cmake --build build/clang-asan --parallel 1
ctest --test-dir build/clang-asan --output-on-failure -E '^cpptest$'
```

On systemd-based Linux, an optional memory scope can protect an interactive workstation; adjust the
limits and parallelism for the machine:

```bash
systemd-run --user --scope \
    -p MemoryHigh=24G -p MemoryMax=32G -p MemorySwapMax=1G \
    cmake --build build/clang-asan --parallel 8
```

The expected validation depth for different changes is listed in
[AGENTS.md](AGENTS.md#validation-cadence).

## Python checks

CTest runs the Python API, adapter, recipe, installation, and mypy tests with the correct build-tree
`PYTHONPATH`. A focused Python test is run like any other CTest test:

```bash
cmake --build build/gcc-debug --target chargefw_python
ctest --test-dir build/gcc-debug --output-on-failure -E '^cpptest$' \
    -R '^test_chargefw_python_calculation$'
```

Ruff is configured in `pyproject.toml`. Run it on modified Python files:

```bash
python -m ruff check path/to/file.py
python -m ruff format --check path/to/file.py
```

## Formatting

Format modified C++ files with the repository configuration:

```bash
clang-format -i --style=file path/to/file.cpp
```

The pre-commit hook formats all tracked C and C++ files:

```bash
pre-commit run --all-files
```

## Installation and containers

Use the self-contained source installation recipes for the [native library](docs/NATIVE.md#build-and-link)
and [CLI](docs/CLI.md#installation). They intentionally exclude the Python and test dependencies used
by the validation presets.

The root `Dockerfile` builds the distributable CLI image:

```bash
docker build --tag chargefw:local .
```

For an image that will remain on the build machine, enable host-native optimization explicitly:

```bash
docker build \
    --build-arg CHARGEFW_ENABLE_NATIVE_OPTIMIZATIONS=ON \
    --tag chargefw:local-native .
```

The Dockerfiles under `docker/` build and test ChargeFW on Ubuntu 26.04, Debian 13, and Fedora 44. Run
these builds one at a time because they are resource intensive.

Python distribution builds and publication are documented in [RELEASING.md](RELEASING.md).
