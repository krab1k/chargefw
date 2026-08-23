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
#include <optional>
#include <stdexcept>
#include <string>
#include <string_view>
#include <vector>

namespace calculation = chargefw::calculation;
namespace charges = chargefw::charges;
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

struct RecordedProgress {
    calculation::CalculationPhase phase{};
    calculation::ExecutionMode mode{};
    std::string method_id;
    std::size_t target_index{};
    std::size_t target_count{};
    std::size_t fragment_index{};
    std::size_t fragment_count{};
    std::size_t molecule_index{};
    std::optional<std::size_t> conformer_index{};
    double elapsed_seconds{};
};

[[nodiscard]] auto snapshot(const calculation::CalculationProgress& progress) -> RecordedProgress {
    return RecordedProgress{.phase = progress.phase,
                            .mode = progress.mode,
                            .method_id = std::string{progress.method_id},
                            .target_index = progress.target_index,
                            .target_count = progress.target_count,
                            .fragment_index = progress.fragment_index,
                            .fragment_count = progress.fragment_count,
                            .molecule_index = progress.molecule_index,
                            .conformer_index = progress.conformer_index,
                            .elapsed_seconds = progress.elapsed_seconds};
}

// Thread-safe recording observer. Captures owned progress snapshots for later assertions.
class RecordingObserver : public calculation::CalculationObserver {
  public:
    void on_progress(const calculation::CalculationProgress& progress) const override {
        const std::scoped_lock lock{mutex_};
        events_.push_back(snapshot(progress));
    }

    [[nodiscard]] auto events() const -> std::vector<RecordedProgress> {
        const std::scoped_lock lock{mutex_};
        return events_;
    }

  private:
    mutable std::mutex mutex_;
    mutable std::vector<RecordedProgress> events_;
};

// Observer that cancels after the first target_started event.
class CancelAfterFirstTarget : public calculation::CalculationObserver {
  public:
    void on_progress(const calculation::CalculationProgress& progress) const override {
        const std::scoped_lock lock{mutex_};
        events_.push_back(snapshot(progress));
        if (progress.phase == calculation::CalculationPhase::target_started) {
            cancel_ = true;
        }
    }

    [[nodiscard]] auto cancelled() const noexcept -> bool override {
        return cancel_;
    }

    [[nodiscard]] auto events() const -> std::vector<RecordedProgress> {
        const std::scoped_lock lock{mutex_};
        return events_;
    }

  private:
    mutable std::mutex mutex_;
    mutable std::vector<RecordedProgress> events_;
    mutable std::atomic<bool> cancel_{false};
};

[[nodiscard]] auto fragment_progress_events(const RecordingObserver& observer)
    -> std::vector<RecordedProgress> {
    auto result = std::vector<RecordedProgress>{};
    for (const auto& event : observer.events()) {
        if (event.phase == calculation::CalculationPhase::fragment_progress) {
            result.push_back(event);
        }
    }
    return result;
}

class ThrowOnTerminalObserver final : public calculation::CalculationObserver {
  public:
    void on_progress(const calculation::CalculationProgress& progress) const override {
        if (progress.phase == calculation::CalculationPhase::computation_finished) {
            terminal_callbacks_.fetch_add(1, std::memory_order_relaxed);
            throw std::runtime_error{"observer failure"};
        }
    }

    [[nodiscard]] auto terminal_callbacks() const noexcept -> std::size_t {
        return terminal_callbacks_.load(std::memory_order_relaxed);
    }

  private:
    mutable std::atomic<std::size_t> terminal_callbacks_{0};
};

auto assert_single_terminal_fragment_progress(const std::vector<RecordedProgress>& events,
                                              const std::size_t fragment_count) -> void {
    assert(!events.empty());
    assert(std::count_if(events.begin(), events.end(), [fragment_count](const auto& event) {
               return event.fragment_index == fragment_count &&
                      event.fragment_count == fragment_count;
           }) == 1);
    for (const auto& event : events) {
        assert(event.fragment_index > 0);
        assert(event.fragment_index <= event.fragment_count);
        assert(event.fragment_count == fragment_count);
    }
}

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

auto make_separated_waters() -> core::Molecule {
    return core::Molecule{std::vector{core::Atom{8}, core::Atom{1}, core::Atom{1}, core::Atom{8},
                                      core::Atom{1}, core::Atom{1}},
                          std::vector{core::Bond{0, 1, core::BondOrder::SINGLE},
                                      core::Bond{0, 2, core::BondOrder::SINGLE},
                                      core::Bond{3, 4, core::BondOrder::SINGLE},
                                      core::Bond{3, 5, core::BondOrder::SINGLE}},
                          {core::Conformer{{{0.0, 0.0, 0.0},
                                            {0.96, 0.0, 0.0},
                                            {-0.24, 0.93, 0.0},
                                            {10.0, 0.0, 0.0},
                                            {10.96, 0.0, 0.0},
                                            {9.76, 0.93, 0.0}}}},
                          "separated-waters"};
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

    // --- Test 8: a terminal observer callback failure does not terminate calculation ---
    {
        const auto observer = ThrowOnTerminalObserver{};
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
        assert(observer.terminal_callbacks() == 1);
    }

    // --- Test 9: serial cutoff execution emits aggregate fragment progress ---
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

        const auto fragment_events = fragment_progress_events(observer);
        assert(!fragment_events.empty());
        assert_single_terminal_fragment_progress(fragment_events,
                                                 fragment_events.front().fragment_count);
    }

    // --- Test 9: serial cover execution emits aggregate multi-pivot progress ---
    {
        const auto observer = RecordingObserver{};
        const auto result = calculate_application(
            calculation::AssessmentRequest{
                .molecules = core::MoleculeCollection{std::vector{make_separated_waters()}},
                .parameter_sets = {},
                .execution_selection =
                    calculation::ExecutionSelection{calculation::ExecutionSelectionKind::cover,
                                                    calculation::minimum_reduced_radius},
                .resource_policy = {.max_threads = 1}},
            observer);

        assert(result.calculated());
        assert(!result.cancelled);

        const auto fragment_events = fragment_progress_events(observer);
        assert(!fragment_events.empty());
        assert_single_terminal_fragment_progress(fragment_events,
                                                 fragment_events.front().fragment_count);
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

    // --- Test 11: empty reduced targets do not emit fragment progress ---
    for (const auto mode : {calculation::ExecutionSelectionKind::cutoff,
                            calculation::ExecutionSelectionKind::cover}) {
        const auto empty_observer = RecordingObserver{};
        const auto empty_molecule = core::Molecule{{}, {}, {core::Conformer{{}}}, "empty"};
        const auto result = calculate_application(
            calculation::AssessmentRequest{
                .molecules = core::MoleculeCollection{std::vector{empty_molecule}},
                .parameter_sets = {},
                .execution_selection =
                    calculation::ExecutionSelection{mode, calculation::minimum_reduced_radius},
                .resource_policy = {.max_threads = 1}},
            empty_observer);
        assert(result.calculated());
        assert(fragment_progress_events(empty_observer).empty());
    }

    // --- Test 12: parallel targets each have one terminal progress snapshot ---
    for (const auto mode : {calculation::ExecutionSelectionKind::cutoff,
                            calculation::ExecutionSelectionKind::cover}) {
        const auto observer = RecordingObserver{};
        const auto result = calculate_application(
            calculation::AssessmentRequest{
                .molecules = core::MoleculeCollection{std::vector{chargefw::test::make_water(),
                                                                  chargefw::test::make_water()}},
                .parameter_sets = {},
                .execution_selection =
                    calculation::ExecutionSelection{mode, calculation::minimum_reduced_radius},
                .resource_policy = {.max_threads = 2}},
            observer);
        assert(result.calculated());

        for (const auto target_index : {std::size_t{0}, std::size_t{1}}) {
            auto target_events = std::vector<RecordedProgress>{};
            for (const auto& event : fragment_progress_events(observer)) {
                if (event.target_index == target_index) {
                    target_events.push_back(event);
                }
            }
            assert(!target_events.empty());
            assert_single_terminal_fragment_progress(target_events,
                                                     target_events.front().fragment_count);
        }
    }

    // --- Test 13: every execution mode preserves source target ordering and identity ---
    for (const auto selection_kind :
         {calculation::ExecutionSelectionKind::full, calculation::ExecutionSelectionKind::cutoff,
          calculation::ExecutionSelectionKind::cover}) {
        const auto observer = RecordingObserver{};
        const auto result = calculate_application(
            calculation::AssessmentRequest{
                .molecules = core::MoleculeCollection{std::vector{
                    chargefw::test::make_two_conformer_water(), chargefw::test::make_water()}},
                .parameter_sets = {},
                .method_id = "eqeq",
                .execution_selection =
                    selection_kind == calculation::ExecutionSelectionKind::full
                        ? calculation::ExecutionSelection{selection_kind}
                        : calculation::ExecutionSelection{selection_kind,
                                                          calculation::minimum_reduced_radius},
                .resource_policy = {.max_threads = 1}},
            observer);
        assert(result.calculated());
        assert(result.charges->size() == 3);

        const auto expected_targets = std::vector<charges::ChargeTarget>{
            {.molecule_index = 0, .conformer_index = 0},
            {.molecule_index = 0, .conformer_index = 1},
            {.molecule_index = 1, .conformer_index = 0},
        };
        for (std::size_t index = 0; index < expected_targets.size(); ++index) {
            assert(result.charges->assignment(index).target.molecule_index ==
                   expected_targets[index].molecule_index);
            assert(result.charges->assignment(index).target.conformer_index ==
                   expected_targets[index].conformer_index);
        }

        auto target_starts = std::vector<RecordedProgress>{};
        for (const auto& event : observer.events()) {
            if (event.phase == calculation::CalculationPhase::target_started) {
                target_starts.push_back(event);
            }
        }
        assert(target_starts.size() == expected_targets.size());
        for (std::size_t index = 0; index < target_starts.size(); ++index) {
            assert(target_starts[index].target_index == index);
            assert(target_starts[index].target_count == expected_targets.size());
            assert(target_starts[index].molecule_index == expected_targets[index].molecule_index);
            assert(target_starts[index].conformer_index == expected_targets[index].conformer_index);
        }
    }

    return 0;
}
