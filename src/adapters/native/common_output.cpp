#include "common_output.h"

#include <chargefw/core/periodic_table.h>

#include <cmath>
#include <format>
#include <fstream>
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

auto bond_type(const core::BondOrder order) -> std::string_view {
    switch (order) {
    case core::BondOrder::SINGLE:
        return "1";
    case core::BondOrder::DOUBLE:
        return "2";
    case core::BondOrder::TRIPLE:
        return "3";
    }
    throw std::invalid_argument{"cannot write unsupported bond order"};
}

auto open_source_file(const std::string& source_path, const std::string_view format_name)
    -> std::ifstream {
    auto input = std::ifstream{source_path, std::ios::binary};
    if (!input) {
        throw std::runtime_error{"unable to open " + std::string{format_name} +
                                 " source file: " + source_path};
    }
    return input;
}

auto validate_assignment(const charges::AtomicCharges& charges, const std::size_t atom_count)
    -> void {
    if (charges.size() != atom_count) {
        throw std::invalid_argument{"charge assignment count does not match molecule atom count"};
    }
}

auto assignment_conformer(const core::Molecule& molecule,
                          const charges::ChargeAssignment& assignment,
                          const std::string_view format_name) -> const core::Conformer& {
    validate_assignment(assignment.charges, molecule.atom_count());
    if (!assignment.target.conformer_index.has_value()) {
        throw std::invalid_argument{"generated " + std::string{format_name} +
                                    " output requires a conformer-specific assignment"};
    }
    const auto conformer_index = *assignment.target.conformer_index;
    if (conformer_index >= molecule.conformer_count()) {
        throw std::invalid_argument{"charge assignment references an unavailable conformer"};
    }
    const auto& conformer = molecule.conformer(conformer_index);
    if (conformer.size() != molecule.atom_count()) {
        throw std::invalid_argument{
            std::string{format_name} +
            " conformer coordinate count does not match molecule atom count"};
    }
    return conformer;
}

} // namespace chargefw::adapters::native::common_output
