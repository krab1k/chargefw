# ChargeFW Delivery Roadmap

This file contains unfinished work only. Implemented state belongs in [PROJECT.md](PROJECT.md), usage
in [README.md](README.md), and implementation rules in [AGENTS.md](AGENTS.md). Remove an item when its
full acceptance criteria are met; do not retain checked history here.

Work is ordered by current product dependency: establish scientific reference evidence, make the native
package consumable and continuously checked, then add Python distribution. ACC III replacement has
separate adoption gates. Optional integrations and distribution formats should be added only for a
concrete supported workflow.

## 1. Scientific validation and ChargeFW2 comparison

The implementation and execution contracts have focused regression coverage. The remaining work is
scientific evidence and release readiness. ChargeFW2 is a previous implementation and useful comparison
material, not an oracle: a difference is a finding to investigate, not evidence by itself that ChargeFW
is wrong. Cutoff and cover remain explicit approximations; the accepted audit establishes their execution
behavior but makes no general reduced-mode accuracy claim.

- [ ] Close the remaining ABEEM parameter-provenance finding in [COMPARISON.md](COMPARISON.md): verify
  the bundled common `k=2.66` against the cited Yang/Shen MEEM source, including its units and derivation.

## 2. Native installation and automation

- [ ] Add a lean CI baseline for the supported development platform: formatting, GCC debug/release
  tests, Clang ASan/UBSan, and the native install/consumer smoke test. Add further compiler, platform,
  or package jobs when they protect a supported distribution rather than maintaining an unused matrix.
- [ ] Give bundled parameter data a machine-readable release identifier and retain it in calculation
  provenance. Before 1.0, document how software, result-schema, and parameter-data versions affect
  compatibility and reproducibility.

## 3. Python and toolkit integration

- [ ] Define a synchronous toolkit-neutral Python API accepting atomic numbers, formal charges,
  indexed bonds, source names/identities, and zero or more coordinate arrays; return NumPy charge
  arrays plus mappings, provenance, and diagnostics.
- [ ] Add nanobind bindings over the owned facade and test ownership, exception translation, array
  validation, collection ordering, source mappings, and all-conformer results.
- [ ] Build a base wheel with bundled parameter data, package-resource discovery that needs no
  environment variable, and clean-install tests for a declared initial CPython/platform set. Expand
  the matrix only to platforms the project intends to support.
- [ ] Add a pure-Python `rdkit.Chem.Mol` converter and an explicit charge-attachment helper behind a lazy
  optional dependency. Preserve atom indices, formal charges, supported bonds, and selected conformers;
  perform no implicit sanitization, hydrogen changes, protonation, embedding, or optimization, and do
  not overwrite existing properties unless requested.

## 4. ACC III adoption gates

- [ ] Establish a reproducible ACC III comparison corpus with licensed/reference inputs, parameter
  versions, options, mappings, and method-specific tolerances. Treat licensed ACC III/ChargeFW2 outputs
  as comparison observations, not automatically correct answers.
- [ ] Investigate each material ChargeFW2/ACC III divergence and classify it as a ChargeFW defect, a
  legacy/reference defect, an intentional change, or unresolved. Retain the rationale and only promote
  behavior to a regression requirement when the scientific or product reason is understood.
- [ ] Before replacing the ACC III backend, define acceptance criteria for its intended workflows and
  obtain independent scientific review of the method comparisons and any cutoff/cover accuracy claims.
