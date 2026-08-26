#include <chargefw/calculation/calculation.h>

#include "calculation/cover_execution.h"
#include "calculation/cutoff_execution.h"
#include "calculation/full_execution.h"
#include "calculation/observer_notifications.h"

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

    ComputationFinishedEmitter(const ComputationFinishedEmitter&) = delete;
    auto operator=(const ComputationFinishedEmitter&) -> ComputationFinishedEmitter& = delete;
    ComputationFinishedEmitter(ComputationFinishedEmitter&&) = delete;
    auto operator=(ComputationFinishedEmitter&&) -> ComputationFinishedEmitter& = delete;

    ~ComputationFinishedEmitter() noexcept {
        detail::report_progress(observer_, CalculationProgress{
                                               .phase = CalculationPhase::computation_finished,
                                               .mode = mode_,
                                               .method_id = method_id_,
                                               .elapsed_seconds =
                                                   std::chrono::duration<double>{
                                                       std::chrono::steady_clock::now() - started_}
                                                       .count(),
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
    if (request.selected.method == nullptr) {
        throw std::invalid_argument{"calculation request has no selected method"};
    }

    const auto computation_started = std::chrono::steady_clock::now();
    const auto computation_finished =
        ComputationFinishedEmitter{request.observer, request.execution_policy.mode(),
                                   request.selected.method->id(), computation_started};
    detail::report_progress(request.observer, CalculationProgress{
                                                  .phase = CalculationPhase::computation_started,
                                                  .mode = request.execution_policy.mode(),
                                                  .method_id = request.selected.method->id(),
                                              });

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
        return ExecutionResult{
            .status = ExecutionStatus::no_executable_plan,
            .charges = std::nullopt,
            .applicability = std::move(assessment.applicability_report_),
            .execution_policy = std::nullopt,
            .execution_issues = {},
            .effective_method_options = std::nullopt,
            .selected_method_id = std::nullopt,
            .selected_parameter_set_id = std::nullopt,
            .failure_message = std::nullopt,
            .metrics = {.applicability_seconds = assessment.applicability_seconds_}};
    }

    if (!assessment.selected_candidate_index_.has_value() ||
        !assessment.execution_policy_.has_value()) {
        throw std::logic_error{"executable calculation assessment has no selected execution plan"};
    }

    const auto& selected =
        assessment.applicability_.applicable.at(*assessment.selected_candidate_index_);
    const auto effective_method_options = selected.method_options;
    const auto selected_method_id = std::string{selected.method->id()};
    const auto selected_parameter_set_id =
        selected.parameter_set == nullptr
            ? std::nullopt
            : std::optional{std::string{selected.parameter_set->id()}};

    const auto computation_started = std::chrono::steady_clock::now();
    auto status = ExecutionStatus::success;
    auto calculated_charges = std::optional<charges::ChargeSet>{};
    auto failure_message = std::optional<std::string>{};

    try {
        auto result =
            calculate(CalculationRequest{.molecules = assessment.prepared_molecules(),
                                         .selected = selected,
                                         .execution_policy = *assessment.execution_policy_,
                                         .max_threads = max_threads,
                                         .observer = observer});
        calculated_charges.emplace(std::move(result.charges));
    } catch (const CalculationCancelled&) {
        status = ExecutionStatus::cancelled;
    } catch (const std::invalid_argument&) {
        throw;
    } catch (const std::exception& error) {
        status = ExecutionStatus::numerical_failure;
        failure_message = error.what();
    }

    const auto computation_seconds =
        std::chrono::duration<double>{std::chrono::steady_clock::now() - computation_started}
            .count();
    return ExecutionResult{.status = status,
                           .charges = std::move(calculated_charges),
                           .applicability = std::move(assessment.applicability_report_),
                           .execution_policy = assessment.execution_policy_,
                           .execution_issues = std::move(assessment.execution_issues_),
                           .effective_method_options = effective_method_options,
                           .selected_method_id = selected_method_id,
                           .selected_parameter_set_id = selected_parameter_set_id,
                           .failure_message = std::move(failure_message),
                           .metrics = {.applicability_seconds = assessment.applicability_seconds_,
                                       .computation_seconds = computation_seconds}};
}

} // namespace chargefw::calculation
