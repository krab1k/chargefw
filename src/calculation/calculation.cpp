#include <chargefw/calculation/calculation.h>

#include "calculation/cover_execution.h"
#include "calculation/cutoff_execution.h"

#include <chargefw/methods/method.h>
#include <chargefw/methods/method_calculation.h>

#include <chrono>
#include <stdexcept>
#include <utility>

namespace chargefw::calculation {

auto calculate(const CalculationRequest& request) -> CalculationResult {
    switch (request.execution_policy.mode()) {
    case ExecutionMode::full:
        return CalculationResult{
            .charges = methods::calculate_charges(request.selected, request.molecules,
                                                  request.max_threads, request.observer)};
    case ExecutionMode::cutoff:
        return CalculationResult{.charges = calculate_cutoff_charges(
                                     request.selected, request.molecules, request.execution_policy,
                                     request.max_threads, request.observer)};
    case ExecutionMode::cover:
        return CalculationResult{.charges = calculate_cover_charges(
                                     request.selected, request.molecules, request.execution_policy,
                                     request.max_threads, request.observer)};
    }

    throw std::invalid_argument{"unknown execution policy"};
}

auto calculate(ApplicationAssessmentResult assessment, const std::size_t max_threads,
               const CalculationObserver* observer) -> ApplicationCalculationResult {
    if (!assessment.executable()) {
        return ApplicationCalculationResult{
            .charges = std::nullopt,
            .applicability = std::move(assessment.applicability),
            .execution_policy = std::nullopt,
            .execution_issues = {},
            .effective_method_options = std::nullopt,
            .metrics = {.applicability_seconds = assessment.applicability_seconds}};
    }

    if (assessment.selected == nullptr || !assessment.execution_policy.has_value()) {
        throw std::logic_error{"executable calculation assessment has no selected execution plan"};
    }

    const auto mode = assessment.execution_policy->mode();
    const auto method_id = assessment.selected->method->id();
    const auto effective_method_options = assessment.selected->method_options;

    if (observer != nullptr) {
        for (const auto& warning : assessment.execution_issues) {
            observer->on_execution_warning(warning);
        }
    }

    const auto computation_started = std::chrono::steady_clock::now();
    if (observer != nullptr) {
        observer->on_progress(CalculationProgress{
            .phase = CalculationPhase::computation_started,
            .mode = mode,
            .method_id = method_id,
        });
    }

    try {
        auto result = calculate(CalculationRequest{.molecules = assessment.prepared_molecules(),
                                                   .selected = *assessment.selected,
                                                   .execution_policy = *assessment.execution_policy,
                                                   .max_threads = max_threads,
                                                   .observer = observer});
        const auto computation_seconds =
            std::chrono::duration<double>{std::chrono::steady_clock::now() - computation_started}
                .count();
        if (observer != nullptr) {
            observer->on_progress(CalculationProgress{
                .phase = CalculationPhase::computation_finished,
                .mode = mode,
                .method_id = method_id,
                .elapsed_seconds = computation_seconds,
            });
        }
        return ApplicationCalculationResult{
            .charges = std::move(result.charges),
            .applicability = std::move(assessment.applicability),
            .execution_policy = assessment.execution_policy,
            .execution_issues = std::move(assessment.execution_issues),
            .effective_method_options = effective_method_options,
            .metrics = {.applicability_seconds = assessment.applicability_seconds,
                        .computation_seconds = computation_seconds}};
    } catch (const CalculationCancelled&) {
        const auto computation_seconds =
            std::chrono::duration<double>{std::chrono::steady_clock::now() - computation_started}
                .count();
        if (observer != nullptr) {
            observer->on_progress(CalculationProgress{
                .phase = CalculationPhase::computation_finished,
                .mode = mode,
                .method_id = method_id,
                .elapsed_seconds = computation_seconds,
            });
        }
        return ApplicationCalculationResult{
            .charges = std::nullopt,
            .applicability = std::move(assessment.applicability),
            .execution_policy = assessment.execution_policy,
            .execution_issues = std::move(assessment.execution_issues),
            .effective_method_options = std::nullopt,
            .metrics = {.applicability_seconds = assessment.applicability_seconds,
                        .computation_seconds = computation_seconds},
            .cancelled = true};
    }
}

auto calculate(const ApplicationCalculationRequest& request) -> ApplicationCalculationResult {
    auto assessment = assess(request);
    if (!assessment.executable() &&
        (request.method_id.has_value() || request.parameter_set_id.has_value() ||
         request.execution_selection.kind() != ExecutionSelectionKind::automatic)) {
        throw std::invalid_argument{"requested calculation selection has no executable plan"};
    }
    return calculate(std::move(assessment), request.resource_policy.max_threads, request.observer);
}

} // namespace chargefw::calculation
