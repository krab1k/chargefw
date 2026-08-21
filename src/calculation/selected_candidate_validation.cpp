#include "calculation/selected_candidate_validation.h"

#include <chargefw/methods/method.h>

#include <stdexcept>
#include <string>

namespace chargefw::calculation::detail {

auto validate_selected_candidate(const methods::ApplicableMethod& selected,
                                 const features::PreparedMoleculeCollection& molecules) -> void {
    if (selected.method == nullptr) {
        throw std::invalid_argument{"selected applicable method has no method"};
    }

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

} // namespace chargefw::calculation::detail
