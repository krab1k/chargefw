#include "support/test_assertions.h"
#include "support/test_molecules.h"
#include "support/test_parameters.h"

#include <algorithm>
#include <atomic>
#include <cassert>
#include <chargefw/calculation/calculation.h>
#include <chargefw/calculation/observer.h>
#include <chargefw/core/molecule.h>
#include <chargefw/core/molecule_collection.h>
#include <chargefw/methods/method_applicability.h>
#include <chargefw/parameters/models/atom_parameters.h>
#include <chargefw/parameters/models/parameter_set.h>
#include <chargefw/parameters/models/parameter_set_metadata.h>
#include <limits>
#include <mutex>
#include <string_view>
#include <vector>

namespace calculation = chargefw::calculation;
namespace core = chargefw::core;
namespace methods = chargefw::methods;

namespace {

[[nodiscard]] auto calculate_application(
    const calculation::AssessmentRequest& request,
    const calculation::CalculationObserver& observer = calculation::default_calculation_observer())
    -> calculation::ExecutionResult {
    const auto max_threads = request.resource_policy.max_threads;
    auto assessment = calculation::assess(request);
    return calculation::calculate(std::move(assessment), max_threads, observer);
}

// Thread-safe recording observer. Captures progress events for later assertions.
class RecordingObserver : public calculation::CalculationObserver {
  public:
    void on_progress(const calculation::CalculationProgress& progress) const override {
        const std::scoped_lock lock{mutex_};
        events_.push_back(progress);
    }

    [[nodiscard]] auto events() const -> std::vector<calculation::CalculationProgress> {
        const std::scoped_lock lock{mutex_};
        return events_;
    }

  private:
    mutable std::mutex mutex_;
    mutable std::vector<calculation::CalculationProgress> events_;
};

// Observer that cancels after the first target_started event.
class CancelAfterFirstTarget : public calculation::CalculationObserver {
  public:
    void on_progress(const calculation::CalculationProgress& progress) const override {
        const std::scoped_lock lock{mutex_};
        events_.push_back(progress);
        if (progress.phase == calculation::CalculationPhase::target_started) {
            cancel_ = true;
        }
    }

    [[nodiscard]] auto cancelled() const noexcept -> bool override {
        return cancel_;
    }

    [[nodiscard]] auto events() const -> std::vector<calculation::CalculationProgress> {
        const std::scoped_lock lock{mutex_};
        return events_;
    }

  private:
    mutable std::mutex mutex_;
    mutable std::vector<calculation::CalculationProgress> events_;
    mutable std::atomic<bool> cancel_{false};
};

auto make_invalid_qeq_parameters() -> chargefw::parameters::ParameterSet {
    return chargefw::parameters::ParameterSet{
        chargefw::parameters::ParameterSetMetadata{
            .id = "invalid-qeq", .method_id = "qeq", .name = "Invalid QEq"},
        {},
        chargefw::parameters::AtomParameters{
            {{.key = chargefw::test::plain_atom_key(1),
              .parameters = {{.name = "electronegativity", .value = 4.5280},
                             {.name = "hardness", .value = 0.0}}},
             {.key = chargefw::test::plain_atom_key(8),
              .parameters = {{.name = "electronegativity", .value = 8.741},
                             {.name = "hardness", .value = 13.364}}}}}};
}

} // namespace

auto main() -> int {
    // --- Test 1: null observer produces the same result as before ---
    {
        const auto result = calculate_application(calculation::AssessmentRequest{
            .molecules = core::MoleculeCollection{std::vector{chargefw::test::make_water()}},
            .parameter_sets = {},
            .method_id = "formal",
            .execution_selection =
                calculation::ExecutionSelection{calculation::ExecutionSelectionKind::full}});

        assert(result.calculated());
        assert(!result.cancelled);
        assert(result.charges->method_id() == std::string_view{"formal"});
    }

    // --- Test 2: computation and target phases are emitted in order ---
    {
        const auto observer = RecordingObserver{};
        const auto result = calculate_application(
            calculation::AssessmentRequest{
                .molecules = core::MoleculeCollection{std::vector{chargefw::test::make_water()}},
                .parameter_sets = {},
                .method_id = "formal",
                .execution_selection =
                    calculation::ExecutionSelection{calculation::ExecutionSelectionKind::full}},
            observer);

        assert(result.calculated());
        assert(!result.cancelled);

        const auto events = observer.events();
        assert(!events.empty());

        // The computation phase starts observation.
        assert(events[0].phase == calculation::CalculationPhase::computation_started);

        // computation_started must carry mode and method_id.
        assert(events[0].mode == calculation::ExecutionMode::full);
        assert(events[0].method_id == std::string_view{"formal"});

        // computation_finished must be the last event and carry mode/method_id.
        assert(events.back().phase == calculation::CalculationPhase::computation_finished);
        assert(events.back().mode == calculation::ExecutionMode::full);
        assert(events.back().method_id == std::string_view{"formal"});
        assert(std::count_if(events.begin(), events.end(), [](const auto& event) {
                   return event.phase == calculation::CalculationPhase::computation_started;
               }) == 1);
        assert(std::count_if(events.begin(), events.end(), [](const auto& event) {
                   return event.phase == calculation::CalculationPhase::computation_finished;
               }) == 1);

        // Between computation_started and computation_finished there must be at least
        // target_started and target_finished.
        const auto computation_start_it =
            std::find_if(events.begin(), events.end(), [](const auto& e) {
                return e.phase == calculation::CalculationPhase::computation_started;
            });
        const auto computation_end_it =
            std::find_if(events.begin(), events.end(), [](const auto& e) {
                return e.phase == calculation::CalculationPhase::computation_finished;
            });

        const auto target_started_it =
            std::find_if(computation_start_it, computation_end_it, [](const auto& e) {
                return e.phase == calculation::CalculationPhase::target_started;
            });
        assert(target_started_it != computation_end_it);

        const auto target_finished_it =
            std::find_if(target_started_it, computation_end_it, [](const auto& e) {
                return e.phase == calculation::CalculationPhase::target_finished;
            });
        assert(target_finished_it != computation_end_it);
    }

    // --- Test 3: assessment does not emit calculation observer events ---
    {
        const auto observer = RecordingObserver{};
        const auto assessment = calculation::assess(calculation::AssessmentRequest{
            .molecules = core::MoleculeCollection{std::vector{chargefw::test::make_water()}},
            .parameter_sets = {},
            .method_id = "formal",
            .execution_selection =
                calculation::ExecutionSelection{calculation::ExecutionSelectionKind::full}});

        assert(assessment.executable());
        assert(observer.events().empty());
    }

    // --- Test 4: threshold warnings remain assessment data before computation ---
    {
        const auto observer = RecordingObserver{};
        auto assessment = calculation::assess(calculation::AssessmentRequest{
            .molecules = core::MoleculeCollection{std::vector{chargefw::test::make_water()}},
            .parameter_sets = {},
            .method_id = "mgc",
            .execution_selection =
                calculation::ExecutionSelection{calculation::ExecutionSelectionKind::full},
            .resource_policy = {.full_atom_threshold = 2}});

        assert(assessment.execution_issues().size() == 1);
        assert(assessment.execution_issues()[0].kind ==
               methods::ExecutionIssueKind::resource_threshold_exceeded);
        assert(observer.events().empty());

        const auto result = calculation::calculate(std::move(assessment), 1, observer);
        assert(result.calculated());
        assert(!result.cancelled);
        assert(observer.events()[0].phase == calculation::CalculationPhase::computation_started);
    }

    // --- Test 5: cancellation produces cancelled result ---
    {
        const auto observer = CancelAfterFirstTarget{};
        const auto result = calculate_application(
            calculation::AssessmentRequest{
                .molecules = core::MoleculeCollection{std::vector{chargefw::test::make_water()}},
                .parameter_sets = {},
                .method_id = "formal",
                .execution_selection =
                    calculation::ExecutionSelection{calculation::ExecutionSelectionKind::full},
                .resource_policy = {.max_threads = 1}},
            observer);

        assert(result.cancelled);
        assert(!result.calculated());
        const auto events = observer.events();
        assert(events.front().phase == calculation::CalculationPhase::computation_started);
        assert(events.back().phase == calculation::CalculationPhase::computation_finished);
        assert(std::count_if(events.begin(), events.end(), [](const auto& event) {
                   return event.phase == calculation::CalculationPhase::computation_finished;
               }) == 1);
    }

    // --- Test 6: validation failures still end observation and propagate unchanged ---
    {
        const auto observer = RecordingObserver{};
        auto assessment = calculation::assess(calculation::AssessmentRequest{
            .molecules = core::MoleculeCollection{std::vector{chargefw::test::make_water()}},
            .parameter_sets = {},
            .method_id = "formal",
            .execution_selection =
                calculation::ExecutionSelection{calculation::ExecutionSelectionKind::full}});

        assert(chargefw::test::throws_invalid_argument([&] -> void {
            static_cast<void>(calculation::calculate(
                std::move(assessment), std::numeric_limits<std::size_t>::max(), observer));
        }));

        const auto events = observer.events();
        assert(events.front().phase == calculation::CalculationPhase::computation_started);
        assert(events.back().phase == calculation::CalculationPhase::computation_finished);
        assert(events.back().mode == calculation::ExecutionMode::full);
        assert(events.back().method_id == std::string_view{"formal"});
        assert(std::count_if(events.begin(), events.end(), [](const auto& event) {
                   return event.phase == calculation::CalculationPhase::computation_finished;
               }) == 1);
    }

    // --- Test 7: solver failures still end observation and propagate unchanged ---
    {
        const auto observer = RecordingObserver{};
        auto assessment = calculation::assess(calculation::AssessmentRequest{
            .molecules = core::MoleculeCollection{std::vector{chargefw::test::make_water()}},
            .parameter_sets = {make_invalid_qeq_parameters()},
            .method_id = "qeq",
            .parameter_set_id = "invalid-qeq",
            .execution_selection =
                calculation::ExecutionSelection{calculation::ExecutionSelectionKind::full}});

        assert(chargefw::test::throws<std::logic_error>([&] -> void {
            static_cast<void>(calculation::calculate(std::move(assessment), 1, observer));
        }));

        const auto events = observer.events();
        assert(events.front().phase == calculation::CalculationPhase::computation_started);
        assert(events.back().phase == calculation::CalculationPhase::computation_finished);
        assert(events.back().mode == calculation::ExecutionMode::full);
        assert(events.back().method_id == std::string_view{"qeq"});
        assert(std::count_if(events.begin(), events.end(), [](const auto& event) {
                   return event.phase == calculation::CalculationPhase::computation_finished;
               }) == 1);
    }

    // --- Test 8: cutoff execution emits fragment progress ---
    {
        const auto observer = RecordingObserver{};
        const auto result = calculate_application(
            calculation::AssessmentRequest{
                .molecules = core::MoleculeCollection{std::vector{chargefw::test::make_water()}},
                .parameter_sets = {},
                .execution_selection =
                    calculation::ExecutionSelection{calculation::ExecutionSelectionKind::cutoff,
                                                    calculation::minimum_reduced_radius},
                .resource_policy = {.max_threads = 1}},
            observer);

        assert(result.calculated());
        assert(!result.cancelled);

        const auto events = observer.events();

        // Must see at least one fragment_finished event.
        const auto fragment_it = std::find_if(events.begin(), events.end(), [](const auto& e) {
            return e.phase == calculation::CalculationPhase::fragment_finished;
        });
        assert(fragment_it != events.end());
        assert(fragment_it->fragment_count > 0);
        // The final fragment event must report full completion.
        const auto last_fragment = std::find_if(events.rbegin(), events.rend(), [](const auto& e) {
            return e.phase == calculation::CalculationPhase::fragment_finished;
        });
        assert(last_fragment->fragment_index == last_fragment->fragment_count);
    }

    // --- Test 9: cover execution emits pivot-selection and fragment progress ---
    {
        const auto observer = RecordingObserver{};
        const auto result = calculate_application(
            calculation::AssessmentRequest{
                .molecules = core::MoleculeCollection{std::vector{chargefw::test::make_water()}},
                .parameter_sets = {},
                .execution_selection =
                    calculation::ExecutionSelection{calculation::ExecutionSelectionKind::cover,
                                                    calculation::minimum_reduced_radius},
                .resource_policy = {.max_threads = 1}},
            observer);

        assert(result.calculated());
        assert(!result.cancelled);

        const auto events = observer.events();

        // Must see at least one fragment_finished event.
        const auto fragment_it = std::find_if(events.begin(), events.end(), [](const auto& e) {
            return e.phase == calculation::CalculationPhase::fragment_finished;
        });
        assert(fragment_it != events.end());
    }

    // --- Test 10: multi-molecule target events carry correct molecule_index ---
    {
        const auto observer = RecordingObserver{};
        const auto result = calculate_application(
            calculation::AssessmentRequest{
                .molecules = core::MoleculeCollection{std::vector{chargefw::test::make_water(),
                                                                  chargefw::test::make_water()}},
                .parameter_sets = {},
                .method_id = "formal",
                .execution_selection =
                    calculation::ExecutionSelection{calculation::ExecutionSelectionKind::full},
                .resource_policy = {.max_threads = 1}},
            observer);

        assert(result.calculated());
        assert(!result.cancelled);

        const auto events = observer.events();

        // Collect all target_started events.
        std::vector<std::size_t> molecule_indices;
        for (const auto& event : events) {
            if (event.phase == calculation::CalculationPhase::target_started) {
                molecule_indices.push_back(event.molecule_index);
            }
        }

        assert(molecule_indices.size() == 2);
        assert(molecule_indices[0] == 0);
        assert(molecule_indices[1] == 1);

        // target_count must be consistent across events.
        for (const auto& event : events) {
            if (event.phase == calculation::CalculationPhase::target_started ||
                event.phase == calculation::CalculationPhase::target_finished) {
                assert(event.target_count == 2);
            }
        }
    }

    return 0;
}
