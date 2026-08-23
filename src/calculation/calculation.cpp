#include <chargefw/calculation/calculation.h>

#include "calculation/cover_execution.h"
#include "calculation/cutoff_execution.h"
#include "calculation/full_execution.h"

#include <chargefw/methods/method.h>

#include <chrono>
#include <stdexcept>
#include <utility>

namespace chargefw::calculation {

namespace {

class ComputationFinishedEmitter {
  public:
    ComputationFinishedEmitter(const CalculationObserver& observer, const ExecutionMode mode,
                               const std::string_view method_id,
                               const std::chrono::steady_clock::time_point started) noexcept
        : observer_{observer}, mode_{mode}, method_id_{method_id}, started_{started} {}

    ~ComputationFinishedEmitter() noexcept {
        observer_.on_progress(CalculationProgress{
            .phase = CalculationPhase::computation_finished,
            .mode = mode_,
            .method_id = method_id_,
            .elapsed_seconds =
                std::chrono::duration<double>{std::chrono::steady_clock::now() - started_}.count(),
        });
    }

  private:
    const CalculationObserver& observer_;
    ExecutionMode mode_;
    std::string_view method_id_;
    std::chrono::steady_clock::time_point started_;
};

} // namespace

auto calculate(const CalculationRequest& request) -> CalculationResult {
    switch (request.execution_policy.mode()) {
    case ExecutionMode::full:
        return CalculationResult{.charges =
                                     calculate_full_charges(request.selected, request.molecules,
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

auto calculate(AssessmentResult assessment, const std::size_t max_threads,
               const CalculationObserver& observer) -> ExecutionResult {
    if (!assessment.executable()) {
        if (assessment.requires_executable_plan()) {
            throw std::invalid_argument{"requested calculation selection has no executable plan"};
        }
        return ExecutionResult{
            .charges = std::nullopt,
            .applicability = std::move(assessment.applicability_report_),
            .execution_policy = std::nullopt,
            .execution_issues = {},
            .effective_method_options = std::nullopt,
            .metrics = {.applicability_seconds = assessment.applicability_seconds_}};
    }

    if (!assessment.selected_candidate_index_.has_value() ||
        !assessment.execution_policy_.has_value()) {
        throw std::logic_error{"executable calculation assessment has no selected execution plan"};
    }

    const auto& selected =
        assessment.applicability_.applicable.at(*assessment.selected_candidate_index_);
    const auto mode = assessment.execution_policy_->mode();
    const auto method_id = selected.method->id();
    const auto effective_method_options = selected.method_options;

    const auto computation_started = std::chrono::steady_clock::now();
    observer.on_progress(CalculationProgress{
        .phase = CalculationPhase::computation_started,
        .mode = mode,
        .method_id = method_id,
    });
    const auto computation_finished =
        ComputationFinishedEmitter{observer, mode, method_id, computation_started};

    try {
        auto result =
            calculate(CalculationRequest{.molecules = assessment.prepared_molecules(),
                                         .selected = selected,
                                         .execution_policy = *assessment.execution_policy_,
                                         .max_threads = max_threads,
                                         .observer = observer});
        const auto computation_seconds =
            std::chrono::duration<double>{std::chrono::steady_clock::now() - computation_started}
                .count();
        return ExecutionResult{
            .charges = std::move(result.charges),
            .applicability = std::move(assessment.applicability_report_),
            .execution_policy = assessment.execution_policy_,
            .execution_issues = std::move(assessment.execution_issues_),
            .effective_method_options = effective_method_options,
            .metrics = {.applicability_seconds = assessment.applicability_seconds_,
                        .computation_seconds = computation_seconds}};
    } catch (const CalculationCancelled&) {
        const auto computation_seconds =
            std::chrono::duration<double>{std::chrono::steady_clock::now() - computation_started}
                .count();
        return ExecutionResult{
            .charges = std::nullopt,
            .applicability = std::move(assessment.applicability_report_),
            .execution_policy = assessment.execution_policy_,
            .execution_issues = std::move(assessment.execution_issues_),
            .effective_method_options = std::nullopt,
            .metrics = {.applicability_seconds = assessment.applicability_seconds_,
                        .computation_seconds = computation_seconds},
            .cancelled = true};
    }
}

} // namespace chargefw::calculation
