#include "methods/applicable_method_execution.h"

#include <chargefw/methods/method.h>

#include <stdexcept>
#include <string>

namespace chargefw::methods::detail {

auto validate_selected_candidate(const ApplicableMethod& selected,
                                 const features::PreparedMoleculeCollection& molecules) -> void {
    if (selected.method == nullptr) {
        throw std::invalid_argument{"selected applicable method has no method"};
    }

    validate_method_options(selected.method->option_schema(), selected.method_options);

    if (selected.method->requires_parameters() != selected.uses_parameters()) {
        throw std::invalid_argument{
            "selected method '" + std::string{selected.method->id()} +
            (selected.method->requires_parameters()
                 ? "' requires parameters, but no parameter set is attached"
                 : "' does not require parameters, but a parameter set is attached")};
    }

    if (!selected.uses_parameters()) {
        if (!selected.classifications.empty()) {
            throw std::invalid_argument{"selected method '" + std::string{selected.method->id()} +
                                        "' has classifications, but no parameter set"};
        }
        return;
    }

    if (selected.classifications.size() != molecules.size()) {
        throw std::invalid_argument{"selected method '" + std::string{selected.method->id()} +
                                    "' has " + std::to_string(selected.classifications.size()) +
                                    " classifications for " + std::to_string(molecules.size()) +
                                    " molecules"};
    }
}

auto validate_coordinate_targets(const ApplicableMethod& selected,
                                 const features::PreparedMoleculeCollection& molecules) -> void {
    if (!selected.method->requirements().coordinates) {
        return;
    }

    for (std::size_t molecule_index = 0; molecule_index < molecules.size(); ++molecule_index) {
        if (molecules[molecule_index].molecule().conformer_count() == 0) {
            throw std::invalid_argument{"selected method '" + std::string{selected.method->id()} +
                                        "' requires coordinates, but molecule " +
                                        std::to_string(molecule_index + 1) + " has no conformers"};
        }
    }
}

auto parameter_set_id_for(const ApplicableMethod& selected) -> std::optional<std::string> {
    return selected.uses_parameters() ? std::optional{std::string{selected.parameter_set->id()}}
                                      : std::nullopt;
}

} // namespace chargefw::methods::detail
