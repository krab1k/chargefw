#include "support/test_assertions.h"
#include "support/test_molecules.h"

#include <algorithm>
#include <atomic>
#include <cassert>
#include <chargefw/calculation/calculation.h>
#include <chargefw/calculation/observer.h>
#include <chargefw/core/molecule.h>
#include <chargefw/core/molecule_collection.h>
#include <chargefw/methods/method_applicability.h>
#include <mutex>
#include <string_view>
#include <vector>

namespace calculation = chargefw::calculation;
namespace core = chargefw::core;
namespace methods = chargefw::methods;

namespace {

// Thread-safe recording observer. Captures progress events and warnings for later assertions.
class RecordingObserver : public calculation::CalculationObserver {
  public:
    void on_progress(const calculation::CalculationProgress& progress) const override {
        const std::scoped_lock lock{mutex_};
        events_.push_back(progress);
        if (progress.phase == calculation::CalculationPhase::computation_started) {
            warnings_precede_computation_ = !warnings_.empty();
        }
    }

    void on_execution_warning(const methods::ExecutionIssue& warning) const override {
        const std::scoped_lock lock{mutex_};
        warnings_.push_back(warning);
    }

    [[nodiscard]] auto events() const -> std::vector<calculation::CalculationProgress> {
        const std::scoped_lock lock{mutex_};
        return events_;
    }

    [[nodiscard]] auto warnings() const -> std::vector<methods::ExecutionIssue> {
        const std::scoped_lock lock{mutex_};
        return warnings_;
    }

    [[nodiscard]] auto warnings_precede_computation() const -> bool {
        const std::scoped_lock lock{mutex_};
        return warnings_precede_computation_;
    }

  private:
    mutable std::mutex mutex_;
    mutable std::vector<calculation::CalculationProgress> events_;
    mutable std::vector<methods::ExecutionIssue> warnings_;
    mutable bool warnings_precede_computation_ = false;
};

// Observer that cancels after the first target_started event.
class CancelAfterFirstTarget : public calculation::CalculationObserver {
  public:
    void on_progress(const calculation::CalculationProgress& progress) const override {
        if (progress.phase == calculation::CalculationPhase::target_started) {
            cancel_ = true;
        }
    }

    [[nodiscard]] auto cancelled() const noexcept -> bool override {
        return cancel_;
    }

  private:
    mutable std::atomic<bool> cancel_{false};
};

} // namespace

auto main() -> int {
    // --- Test 1: null observer produces the same result as before ---
    {
        const auto result = calculation::calculate(calculation::ApplicationCalculationRequest{
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
        const auto result = calculation::calculate(calculation::ApplicationCalculationRequest{
            .molecules = core::MoleculeCollection{std::vector{chargefw::test::make_water()}},
            .parameter_sets = {},
            .method_id = "formal",
            .execution_selection =
                calculation::ExecutionSelection{calculation::ExecutionSelectionKind::full},
            .observer = &observer});

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
        const auto assessment = calculation::assess(calculation::ApplicationCalculationRequest{
            .molecules = core::MoleculeCollection{std::vector{chargefw::test::make_water()}},
            .parameter_sets = {},
            .method_id = "formal",
            .execution_selection =
                calculation::ExecutionSelection{calculation::ExecutionSelectionKind::full},
            .observer = &observer});

        assert(assessment.executable());
        assert(observer.events().empty());
        assert(observer.warnings().empty());
    }

    // --- Test 4: threshold warnings are reported before computation ---
    {
        const auto observer = RecordingObserver{};
        const auto result = calculation::calculate(calculation::ApplicationCalculationRequest{
            .molecules = core::MoleculeCollection{std::vector{chargefw::test::make_water()}},
            .parameter_sets = {},
            .method_id = "mgc",
            .execution_selection =
                calculation::ExecutionSelection{calculation::ExecutionSelectionKind::full},
            .resource_policy = {.full_atom_threshold = 2},
            .observer = &observer});

        assert(result.calculated());
        assert(!result.cancelled);
        assert(result.execution_issues.size() == 1);
        assert(result.execution_issues[0].kind ==
               methods::ExecutionIssueKind::resource_threshold_exceeded);

        const auto warnings = observer.warnings();
        assert(warnings.size() == 1);
        assert(warnings[0].kind == methods::ExecutionIssueKind::resource_threshold_exceeded);
        assert(observer.warnings_precede_computation());
    }

    // --- Test 5: cancellation produces cancelled result ---
    {
        const auto observer = CancelAfterFirstTarget{};
        const auto result = calculation::calculate(calculation::ApplicationCalculationRequest{
            .molecules = core::MoleculeCollection{std::vector{chargefw::test::make_water()}},
            .parameter_sets = {},
            .method_id = "formal",
            .execution_selection =
                calculation::ExecutionSelection{calculation::ExecutionSelectionKind::full},
            .resource_policy = {.max_threads = 1},
            .observer = &observer});

        assert(result.cancelled);
        assert(!result.calculated());
    }

    // --- Test 6: cutoff execution emits fragment progress ---
    {
        const auto observer = RecordingObserver{};
        const auto result = calculation::calculate(calculation::ApplicationCalculationRequest{
            .molecules = core::MoleculeCollection{std::vector{chargefw::test::make_water()}},
            .parameter_sets = {},
            .execution_selection =
                calculation::ExecutionSelection{calculation::ExecutionSelectionKind::cutoff,
                                                calculation::minimum_reduced_radius},
            .resource_policy = {.max_threads = 1},
            .observer = &observer});

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

    // --- Test 7: cover execution emits pivot-selection and fragment progress ---
    {
        const auto observer = RecordingObserver{};
        const auto result = calculation::calculate(calculation::ApplicationCalculationRequest{
            .molecules = core::MoleculeCollection{std::vector{chargefw::test::make_water()}},
            .parameter_sets = {},
            .execution_selection =
                calculation::ExecutionSelection{calculation::ExecutionSelectionKind::cover,
                                                calculation::minimum_reduced_radius},
            .resource_policy = {.max_threads = 1},
            .observer = &observer});

        assert(result.calculated());
        assert(!result.cancelled);

        const auto events = observer.events();

        // Must see at least one fragment_finished event.
        const auto fragment_it = std::find_if(events.begin(), events.end(), [](const auto& e) {
            return e.phase == calculation::CalculationPhase::fragment_finished;
        });
        assert(fragment_it != events.end());
    }

    // --- Test 8: multi-molecule target events carry correct molecule_index ---
    {
        const auto observer = RecordingObserver{};
        const auto result = calculation::calculate(calculation::ApplicationCalculationRequest{
            .molecules = core::MoleculeCollection{std::vector{chargefw::test::make_water(),
                                                              chargefw::test::make_water()}},
            .parameter_sets = {},
            .method_id = "formal",
            .execution_selection =
                calculation::ExecutionSelection{calculation::ExecutionSelectionKind::full},
            .resource_policy = {.max_threads = 1},
            .observer = &observer});

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
