#include "common_output.h"

#include <chargefw/core/periodic_table.h>

#include <cmath>
#include <format>
#include <stdexcept>
#include <string>

namespace chargefw::adapters::native::common_output {
namespace {

constexpr auto charge_scale = 10000.0;

} // namespace

auto formatted_charge(const double value) -> std::string {
    return std::format("{:.4f}", std::round(value * charge_scale) / charge_scale);
}

auto generated_atom_name(const core::Atom& atom, const std::size_t index) -> std::string {
    if (!atom.name().empty()) {
        return std::string{atom.name()};
    }
    return std::string{core::element_symbol(atom.atomic_number())} + std::to_string(index + 1);
}

auto atom_element_symbol(const core::Atom& atom) -> std::string {
    return std::string{core::element_symbol(atom.atomic_number())};
}

auto bond_type(const core::BondOrder order, const ::chargefw::adapters::native::BondFormat format)
    -> std::string_view {
    switch (order) {
    case core::BondOrder::SINGLE:
        return "1";
    case core::BondOrder::DOUBLE:
        return "2";
    case core::BondOrder::TRIPLE:
        return "3";
    case core::BondOrder::AROMATIC:
        return format == ::chargefw::adapters::native::BondFormat::mol ? "4" : "ar";
    case core::BondOrder::UNKNOWN:
        break;
    }
    throw std::invalid_argument{format == ::chargefw::adapters::native::BondFormat::mol
                                    ? "cannot write MOL/SDF unknown bond order"
                                    : "cannot write MOL2 unknown bond order"};
}

auto validate_assignment(const charges::AtomicCharges& charges, const std::size_t atom_count)
    -> void {
    if (charges.size() != atom_count) {
        throw std::invalid_argument{"charge assignment count does not match molecule atom count"};
    }
}

} // namespace chargefw::adapters::native::common_output
