#pragma once

#include <chargefw/core/molecule.h>

#include <cstddef>
#include <span>

namespace chargefw::test {

/// Water graph (OH2) without coordinates. Use for topology-only methods, where conformers
/// are semantically irrelevant.
[[nodiscard]] auto make_water_graph() -> core::Molecule;

/// Water with a single conformer. Use for geometry-dependent methods.
[[nodiscard]] auto make_water() -> core::Molecule;

/// Water with two distinct conformers. Use to verify per-conformer fan-out of
/// geometry-dependent methods.
[[nodiscard]] auto make_two_conformer_water() -> core::Molecule;

[[nodiscard]] auto make_formally_charged_pair() -> core::Molecule;

/// Returns a numbering-only relabeling of the molecule: atom `new_atom_order[new_index]` of the
/// original becomes atom `new_index`. Bonds are remapped accordingly and conformer positions
/// follow the same permutation, so the physical molecule is unchanged. For a deterministic
/// method, charge `i` of the relabeled result must equal charge `new_atom_order[i]` of the
/// original result (up to floating-point associativity).
[[nodiscard]] auto relabel_atoms(const core::Molecule& molecule,
                                 std::span<const std::size_t> new_atom_order) -> core::Molecule;

/// Returns the molecule with every bond endpoint swapped. Bond endpoints carry no direction,
/// so any method must produce identical charges.
[[nodiscard]] auto flip_bond_directions(const core::Molecule& molecule) -> core::Molecule;

} // namespace chargefw::test
