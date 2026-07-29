#include "common_output.h"

#include <chargefw/core/periodic_table.h>

#include <cmath>
#include <iomanip>
#include <sstream>
#include <stdexcept>
#include <string>

namespace chargefw::adapters::native::common_output {
namespace {

constexpr auto charge_scale = 10000.0;

} // namespace

auto formatted_charge(const double value) -> std::string {
    auto output = std::ostringstream{};
    output << std::fixed << std::setprecision(4) << std::round(value * charge_scale) / charge_scale;
    return output.str();
}

auto generated_atom_name(const core::Atom& atom, const std::size_t index) -> std::string {
    if (!atom.name().empty()) {
        return std::string{atom.name()};
    }
    return std::string{core::element_symbol(atom.atomic_number())} + std::to_string(index + 1);
}

auto mol2_atom_type(const core::Atom& atom) -> std::string {
    return std::string{core::element_symbol(atom.atomic_number())};
}

auto mol2_bond_type(const core::BondOrder order) -> std::string {
    switch (order) {
    case core::BondOrder::SINGLE:
        return "1";
    case core::BondOrder::DOUBLE:
        return "2";
    case core::BondOrder::TRIPLE:
        return "3";
    case core::BondOrder::AROMATIC:
        return "ar";
    case core::BondOrder::UNKNOWN:
        break;
    }
    throw std::invalid_argument{"cannot write MOL2 unknown bond order"};
}

auto validate_assignment(const charges::AtomicCharges& charges, const std::size_t atom_count)
    -> void {
    if (charges.size() != atom_count) {
        throw std::invalid_argument{"charge assignment count does not match molecule atom count"};
    }
}

} // namespace chargefw::adapters::native::common_output
