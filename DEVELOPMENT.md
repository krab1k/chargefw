# ChargeFW development

This is the practical guide to building and testing ChargeFW. Repository policy and architecture rules
are in [AGENTS.md](AGENTS.md); public behavior is documented under [`docs/`](docs/).

## Requirements

Native development requires CMake 3.28 or newer, Ninja, and a C++23-capable GCC or Clang toolchain.
CMake downloads pinned dependencies when compatible system packages are unavailable.

The checked-in presets also build the Python bindings. Their Python interpreter needs Python 3.10 or
newer, development headers, NumPy, Nanobind, mypy, and Gemmi 0.7.4. RDKit is optional. Override the
preset's `/usr/bin/python3` when these packages are installed elsewhere:

```bash
cmake --preset gcc-debug \
    -DCHARGEFW_PYTHON_EXECUTABLE=/path/to/python
```

## Build and test

Use `gcc-debug` for normal development:

```bash
cmake --preset gcc-debug
cmake --build --preset gcc-debug
ctest --preset gcc-debug
```

For a quick edit cycle, build the affected target and run its test by name:

```bash
cmake --build --preset gcc-debug --target test_molecule
ctest --preset gcc-debug -R '^test_molecule$'
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

Configure, build, and test another preset in the same way. `clang-tidy` runs during its build and has no
CTest preset:

```bash
cmake --preset clang-tidy
cmake --build --preset clang-tidy
```

Build sanitizer presets with low parallelism because they can use substantial memory:

```bash
cmake --preset clang-asan
cmake --build --preset clang-asan -j 1
ctest --preset clang-asan
```

The expected validation depth for different changes is listed in
[AGENTS.md](AGENTS.md#validation-cadence).

## Python checks

CTest runs the Python API, adapter, recipe, installation, and mypy tests with the correct build-tree
`PYTHONPATH`. A focused Python test is run like any other CTest test:

```bash
cmake --build --preset gcc-debug --target chargefw_python
ctest --preset gcc-debug -R '^test_chargefw_python_calculation$'
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

Build and install the native library and CLI:

```bash
cmake --preset gcc-release -DCMAKE_INSTALL_PREFIX="$PWD/_install"
cmake --build --preset gcc-release
cmake --install build/gcc-release --strip
```

The root `Dockerfile` builds the distributable CLI image:

```bash
docker build --tag chargefw:local .
```

The Dockerfiles under `docker/` build and test ChargeFW on Ubuntu 26.04, Debian 13, and Fedora 44. Run
these builds one at a time because they are resource intensive.

Python distribution builds and publication are documented in [RELEASING.md](RELEASING.md).
