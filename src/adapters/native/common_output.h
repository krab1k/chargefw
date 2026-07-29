#pragma once

#include <chargefw/charges/atomic_charges.h>
#include <chargefw/core/atom.h>
#include <chargefw/core/bond.h>

#include <cstddef>
#include <string>

namespace chargefw::adapters::native::common_output {

[[nodiscard]] auto formatted_charge(double value) -> std::string;
[[nodiscard]] auto generated_atom_name(const core::Atom& atom, std::size_t index) -> std::string;
[[nodiscard]] auto mol2_atom_type(const core::Atom& atom) -> std::string;
[[nodiscard]] auto mol2_bond_type(core::BondOrder order) -> std::string;

auto validate_assignment(const charges::AtomicCharges& charges, std::size_t atom_count) -> void;

} // namespace chargefw::adapters::native::common_output
