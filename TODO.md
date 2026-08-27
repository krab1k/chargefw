# ChargeFW Delivery Roadmap

This file contains unfinished work only. Implemented state belongs in [PROJECT.md](PROJECT.md), usage
in [README.md](README.md), and implementation rules in [AGENTS.md](AGENTS.md). Remove an item when its
full acceptance criteria are met; do not retain checked history here.

Work is ordered by current product dependency: add Python distribution, map the resulting API to the
ACC III backend's actual needs, then automate release validation. Optional integrations and distribution
formats should be added only for a concrete supported workflow.

## 1. Python and toolkit integration

The detailed API contract, packaging decisions, adapter boundaries, and implementation sequence are
maintained in [PYTHON.md](PYTHON.md); this roadmap retains only release-level acceptance criteria.

- [ ] Define a synchronous toolkit-neutral Python API accepting atomic numbers, formal charges,
  indexed bonds, source names/identities, and zero or more coordinate arrays; return NumPy charge
  arrays plus mappings, provenance, and diagnostics.
- [ ] Add nanobind bindings over the owned facade and test ownership, exception translation, array
  validation, collection ordering, source mappings, and all-conformer results.
- [ ] Build a base wheel with bundled parameter data, package-resource discovery that needs no
  environment variable, required tested Gemmi Python integration, and clean-install tests for a declared
  initial CPython/platform set. Expand the matrix only to platforms the project intends to support.
- [ ] Add a pure-Python `rdkit.Chem.Mol` converter and an explicit charge-attachment helper behind a lazy
  optional dependency. Preserve atom indices, formal charges, supported bonds, and selected conformers;
  perform no implicit sanitization, hydrogen changes, protonation, embedding, or optimization, and do
  not overwrite existing properties unless requested.

## 2. ACC III Python-backend integration

- [ ] Map the ChargeFW2 capabilities actually invoked by ACC III—methods, parameter sets, options,
  molecule and conformer inputs, source mappings, outputs, and failures—and provide the corresponding
  explicit capabilities through the ChargeFW Python API/bindings. Document any intentionally unsupported
  legacy behavior before migrating the ACC III backend.

## 3. Release hardening and automation

- [ ] Before a native/package release, automate the existing native validation baseline in CI for release
  candidates and protected maintenance branches: formatting verification, GCC debug/release CTest suites,
  Clang ASan/UBSan CTest suites, and the installed-package/downstream-consumer relocation smoke tests.
  The presets and smoke tests already exist and remain the normal local development workflow. Add further
  compiler, platform, or package jobs only when they protect a supported distribution.
