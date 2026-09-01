# ChargeFW TODO

This file lists unfinished product work. Implemented behavior is documented under [`docs/`](docs/).

## Distribution

- [ ] Declare and qualify an initial CPython/Linux wheel matrix. Each wheel must pass clean-environment
  installation, relocation, bundled-parameter discovery, shared-library loading, Gemmi conversion, and
  calculation tests without relying on another ChargeFW installation or `LD_LIBRARY_PATH`.
- [ ] Automate release validation in CI, including formatting, GCC and Clang debug/release tests,
  AddressSanitizer and UndefinedBehaviorSanitizer tests, and installed-package relocation and downstream
  CMake consumer tests.

## Integrations

- [ ] Add a lazy, optional RDKit adapter that converts `rdkit.Chem.Mol` values without implicit
  sanitization, hydrogen changes, protonation, embedding, or optimization. Preserve atom, bond, formal
  charge, and conformer mappings, and provide an explicit charge-attachment helper with overwrite
  protection.
- [ ] Compare the Python API with the capabilities actually used by the ACC III backend and add only the
  missing explicit capabilities needed for migration. Document intentionally unsupported legacy
  behavior before replacing the existing backend.

## Scientific qualification

- [ ] Establish method-specific validation and documented tolerances before making broad accuracy claims
  for cutoff or cover execution. Current whole-molecule-radius and representative-structure checks verify
  implementation behavior, not a general reduced-mode accuracy envelope.
