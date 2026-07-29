#include <chargefw/calculation/calculation.h>

#include <chargefw/methods/method.h>
#include <chargefw/methods/method_calculation.h>
#include <chargefw/methods/method_registry.h>

#include <algorithm>
#include <stdexcept>
#include <string>
#include <string_view>
#include <vector>

namespace chargefw::calculation {
namespace {

[[nodiscard]] auto parameter_priority_of(const methods::ApplicableMethod& candidate) noexcept
    -> unsigned int {
    if (candidate.parameter_set == nullptr) {
        return 0;
    }

    return candidate.parameter_set->priority();
}

[[nodiscard]] auto parameter_id_of(const methods::ApplicableMethod& candidate) noexcept
    -> std::string_view {
    if (candidate.parameter_set == nullptr) {
        return {};
    }

    return candidate.parameter_set->id();
}

[[nodiscard]] auto ranks_before(const methods::ApplicableMethod& first,
                                const methods::ApplicableMethod& second) noexcept -> bool {
    const auto first_method_priority = first.method->metadata().priority;
    const auto second_method_priority = second.method->metadata().priority;

    if (first_method_priority != second_method_priority) {
        return first_method_priority > second_method_priority;
    }

    const auto first_parameter_priority = parameter_priority_of(first);
    const auto second_parameter_priority = parameter_priority_of(second);

    if (first_parameter_priority != second_parameter_priority) {
        return first_parameter_priority > second_parameter_priority;
    }

    if (first.method->id() != second.method->id()) {
        return first.method->id() < second.method->id();
    }

    return parameter_id_of(first) < parameter_id_of(second);
}

[[nodiscard]] auto application_methods(const ApplicationCalculationRequest& request)
    -> std::vector<const methods::Method*> {
    const auto& registry = methods::method_registry();

    if (request.method_id.has_value()) {
        const auto* method = registry.find(*request.method_id);

        if (method == nullptr) {
            throw std::invalid_argument{"method '" + *request.method_id + "' is not registered"};
        }

        return {method};
    }

    std::vector<const methods::Method*> result;
    result.reserve(registry.methods().size());

    for (const auto& method : registry.methods()) {
        result.push_back(method.get());
    }

    return result;
}

[[nodiscard]] auto application_parameter_sets(const ApplicationCalculationRequest& request)
    -> std::vector<parameters::ParameterSet> {
    if (!request.parameter_set_id.has_value()) {
        return request.parameter_sets;
    }

    const auto found = std::ranges::find_if(
        request.parameter_sets, [&request](const parameters::ParameterSet& parameter_set) -> bool {
            return parameter_set.id() == *request.parameter_set_id;
        });

    if (found == request.parameter_sets.end()) {
        throw std::invalid_argument{"parameter set '" + *request.parameter_set_id +
                                    "' was not provided"};
    }

    return {*found};
}

} // namespace

auto calculate(const CalculationRequest& request) -> CalculationResult {
    auto applicability = methods::find_applicable_methods(
        request.molecules, request.candidate_methods, request.parameter_sets);

    if (applicability.empty()) {
        return CalculationResult{.charges = std::nullopt,
                                 .applicability = std::move(applicability)};
    }

    const auto selected = std::ranges::min_element(applicability.applicable, ranks_before);

    return CalculationResult{.charges = methods::calculate_charges(*selected, request.molecules),
                             .applicability = std::move(applicability)};
}

auto calculate(const ApplicationCalculationRequest& request) -> ApplicationCalculationResult {
    const auto candidate_methods = application_methods(request);
    const auto parameter_sets = application_parameter_sets(request);
    const features::PreparedMoleculeCollection prepared{request.molecules};

    auto result = calculate(CalculationRequest{.molecules = prepared,
                                               .candidate_methods = candidate_methods,
                                               .parameter_sets = parameter_sets});

    if (!result.calculated() &&
        (request.method_id.has_value() || request.parameter_set_id.has_value())) {
        throw std::invalid_argument{"requested calculation selection is not applicable"};
    }

    return result;
}

} // namespace chargefw::calculation
