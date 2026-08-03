#pragma once

#include <chargefw/adapters/gemmi/input_options.h>

#include <gemmi/model.hpp>

#include <cstddef>

namespace chargefw::adapters::gemmi::selection {

[[nodiscard]] auto include_residue(const ::gemmi::Residue& residue, RecordSelection selection)
    -> bool;

[[nodiscard]] auto is_first_named_atom(const ::gemmi::Residue& residue, std::size_t atom_index)
    -> bool;

[[nodiscard]] auto select_altloc(const ::gemmi::Residue& residue, std::size_t first_atom_index)
    -> const ::gemmi::Atom&;

} // namespace chargefw::adapters::gemmi::selection
