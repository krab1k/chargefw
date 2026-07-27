#include <chargefw/calculation/calculation.h>

#include <chargefw/methods/method.h>
#include <chargefw/methods/method_calculation.h>

#include <algorithm>
#include <string_view>

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

} // namespace chargefw::calculation
