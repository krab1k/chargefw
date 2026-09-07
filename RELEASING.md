# ChargeFW releases

This is the manual procedure for building and publishing the Python package. Development builds and the
full native validation matrix are covered by [DEVELOPMENT.md](DEVELOPMENT.md).

## Current wheel matrix

The current manual release target is CPython 3.10 through 3.14 on Linux x86-64. `cibuildwheel` builds each
wheel in a manylinux container, repairs its native libraries, installs it, and runs the Python test suite.
The resulting wheels carry `manylinux_2_27_x86_64` and `manylinux_2_28_x86_64` tags.

CPython 3.15 is not included because the pinned optional Gemmi 0.7.4 package does not currently provide a
usable CPython 3.15 wheel. Other operating systems and architectures are not yet qualified.

## Before building

The release version comes from the top-level `project(... VERSION ...)` declaration in
`CMakeLists.txt`. It is also used for the native library and `chargefw.__version__`.

Before a release:

1. Update the version in `CMakeLists.txt`.
2. Run the required debug, release, sanitizer, and static-analysis checks.
3. Commit the release changes and start from a clean worktree.
4. Install `uv` and either Podman or Docker.

Published files cannot be replaced. Use a new version whenever artifacts must be rebuilt after upload.

## Build

Build the source distribution first, then build every wheel from that exact archive:

```bash
VERSION=0.1.1  # Must match CMakeLists.txt.

uv build --sdist --out-dir release-wheelhouse --clear

CIBW_CONTAINER_ENGINE=podman uvx cibuildwheel==3.4.1 \
    --platform linux \
    --output-dir release-wheelhouse \
    "release-wheelhouse/chargefw-$VERSION.tar.gz"
```

Omit `CIBW_CONTAINER_ENGINE=podman` to use Docker. The wheel matrix, tests, and portable CMake settings
are defined in `pyproject.toml`.

A successful build produces one source distribution and five wheels in `release-wheelhouse/`.

## Check artifacts

Validate all package metadata and record checksums:

```bash
uvx twine check --strict release-wheelhouse/*
sha256sum release-wheelhouse/* > /tmp/chargefw-"$VERSION".sha256
```

`cibuildwheel` has already installed and tested every wheel in an isolated container. Before uploading,
also inspect the filenames, version, Python tags, and manylinux tags. Keep the checksum file until the
release is complete.

## Publish to TestPyPI

Set `TWINE_USERNAME=__token__` and `TWINE_PASSWORD` to a TestPyPI project token. Keep credentials outside
the repository and shell history:

```bash
uvx twine upload --repository testpypi release-wheelhouse/*
```

Install the uploaded package in a clean environment and run a smoke calculation:

```bash
uv venv --clear /tmp/chargefw-testpypi --python 3.14
uv pip install --python /tmp/chargefw-testpypi/bin/python 'numpy>=1.26'
uv pip install \
    --python /tmp/chargefw-testpypi/bin/python \
    --no-deps \
    --no-cache \
    --only-binary chargefw \
    --default-index https://test.pypi.org/simple/ \
    "chargefw==$VERSION"

/tmp/chargefw-testpypi/bin/python docs/recipes/calculate_file.py \
    tests/fixtures/synthetic/sdf/water.sdf \
    --format sdf \
    --method qeq \
    --parameter-set QEq_original
```

## Publish to PyPI

Publish the same files that passed TestPyPI validation; do not rebuild them. Set the Twine variables to a
production PyPI project token. Compare the TestPyPI file hashes with
`/tmp/chargefw-$VERSION.sha256`; stop if any file differs.

Create the tag locally, then upload:

```bash
git tag -a "v$VERSION" -m "ChargeFW $VERSION"
uvx twine upload release-wheelhouse/*
```

After publication, install the exact version from PyPI in a clean environment, repeat the smoke test, and
push the tag:

```bash
git push origin "v$VERSION"
```

If the original TestPyPI-approved files or checksums are no longer available, do not substitute rebuilt
files under the same version. Increment the version and repeat the release process.

Release automation and additional wheel platforms are tracked in [TODO.md](TODO.md#distribution).
