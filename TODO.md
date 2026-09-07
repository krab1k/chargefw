# ChargeFW TODO

This file lists unfinished product work. Implemented behavior is documented under [`docs/`](docs/), and
developer and release procedures are in [DEVELOPMENT.md](DEVELOPMENT.md) and
[RELEASING.md](RELEASING.md).

## Distribution

- [ ] Automate the manual `cibuildwheel` release process in [RELEASING.md](RELEASING.md), including
  artifact retention and publication through trusted PyPI environments.
- [ ] Automate repository validation, including formatting, GCC and Clang debug/release tests,
  AddressSanitizer and UndefinedBehaviorSanitizer tests, installed-package relocation, and the downstream
  CMake consumer.
- [ ] Add and qualify wheels for other platforms, architectures, and Python versions as their native and
  optional Python dependencies permit.

## Integrations

- [ ] Qualify the optional RDKit adapter against supported RDKit releases with real-toolkit conversion,
  mapping, property attachment, and SD serialization tests.
- [ ] Compare the Python API with the capabilities actually used by the ACC III backend and add only the
  missing explicit capabilities needed for migration. Document intentionally unsupported legacy
  behavior before replacing the existing backend.

## Scientific qualification

- [ ] Establish method-specific validation and documented tolerances before making broad accuracy claims
  for cutoff or cover execution. Current whole-molecule-radius and representative-structure checks verify
  implementation behavior, not a general reduced-mode accuracy envelope.
