#pragma once

#include <chargefw/charges/atomic_charges.h>
#include <chargefw/charges/charge_collection.h>
#include <chargefw/core/atom.h>
#include <chargefw/core/bond.h>
#include <chargefw/core/molecule.h>

#include "bond_format.h"

#include <cstddef>
#include <fstream>
#include <string>
#include <string_view>

namespace chargefw::adapters::native::common_output {

[[nodiscard]] auto formatted_charge(double value) -> std::string;
[[nodiscard]] auto generated_atom_name(const core::Atom& atom, std::size_t index) -> std::string;
[[nodiscard]] auto atom_element_symbol(const core::Atom& atom) -> std::string;
[[nodiscard]] auto bond_type(core::BondOrder order, ::chargefw::adapters::native::BondFormat format)
    -> std::string_view;

[[nodiscard]] auto open_source_file(const std::string& source_path, std::string_view format_name)
    -> std::ifstream;

auto validate_assignment(const charges::AtomicCharges& charges, std::size_t atom_count) -> void;

[[nodiscard]] auto assignment_conformer(const core::Molecule& molecule,
                                        const charges::ChargeAssignment& assignment,
                                        std::string_view format_name) -> const core::Conformer&;

} // namespace chargefw::adapters::native::common_output
