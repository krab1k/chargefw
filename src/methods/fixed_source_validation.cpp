#include "fixed_source_validation.h"

#include <cmath>
#include <stdexcept>
#include <string>

namespace chargefw::methods::detail {

auto validate_fixed_source_input(const CalculationInput& input, const std::string_view method_name)
    -> void {
    const auto sources = input.fixed_sources();
    if (sources.empty()) {
        return;
    }

    const auto method = std::string{method_name};
    if (!std::isfinite(input.target_charge())) {
        throw std::invalid_argument{method +
                                    " fixed-source calculation requires a finite target charge"};
    }
    const auto& molecule = input.molecule();
    if (molecule.atom_count() == 0) {
        throw std::invalid_argument{method + " fixed-source calculation requires active atoms"};
    }
    for (std::size_t source_index = 0; source_index < sources.size(); ++source_index) {
        const auto& source = sources[source_index];
        if (!std::isfinite(source.charge) || !std::isfinite(source.position.x) ||
            !std::isfinite(source.position.y) || !std::isfinite(source.position.z)) {
            throw std::invalid_argument{method + " fixed source " + std::to_string(source_index) +
                                        " has non-finite charge or coordinates"};
        }
    }

    const auto& geometry = input.geometry();
    if (&geometry.molecule() != &molecule) {
        throw std::invalid_argument{
            method + " fixed-source geometry does not belong to the active molecule"};
    }
    if (const auto atom_index = geometry.first_nonfinite_atom_index(); atom_index.has_value()) {
        throw std::invalid_argument{method +
                                    " fixed-source calculation has non-finite coordinates at "
                                    "active atom " +
                                    std::to_string(*atom_index)};
    }
}

} // namespace chargefw::methods::detail
