#include "support/test_molecules.h"
#include "support/test_parameters.h"

#include <algorithm>
#include <atomic>
#include <chargefw/calculation/calculation.h>
#include <chargefw/calculation/observer.h>
#include <chargefw/core/molecule.h>
#include <chargefw/core/molecule_collection.h>
#include <chargefw/features/prepared_molecule_collection.h>
#include <chargefw/methods/method.h>
#include <chargefw/methods/method_applicability.h>
#include <chargefw/methods/method_metadata.h>
#include <chargefw/parameters/models/atom_parameters.h>
#include <chargefw/parameters/models/parameter_set.h>
#include <chargefw/parameters/models/parameter_set_metadata.h>
#include <limits>
#include <mutex>
#include <optional>
#include <snitch/snitch.hpp>
#include <stdexcept>
#include <string>
#include <string_view>
#include <vector>

namespace calculation = chargefw::calculation;
namespace charges = chargefw::charges;
namespace core = chargefw::core;
namespace features = chargefw::features;
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
    std::size_t completed_fragment_count{};
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
                            .completed_fragment_count = progress.completed_fragment_count,
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

// Observer that waits for completed reduced fragment work before requesting cancellation. The
// post-request cancelled() counter proves cancellation was observed by a subsequent fragment-loop
// iteration rather than by the target-boundary check.
class CancelAfterFirstFragmentProgress final : public calculation::CalculationObserver {
  public:
    void on_progress(const calculation::CalculationProgress& progress) const override {
        const std::scoped_lock lock{mutex_};
        events_.push_back(snapshot(progress));
        if (progress.phase == calculation::CalculationPhase::fragment_progress) {
            cancel_.store(true, std::memory_order_relaxed);
        }
    }

    [[nodiscard]] auto cancelled() const noexcept -> bool override {
        if (cancel_.load(std::memory_order_relaxed)) {
            cancellation_checks_after_request_.fetch_add(1, std::memory_order_relaxed);
            return true;
        }
        return false;
    }

    [[nodiscard]] auto events() const -> std::vector<RecordedProgress> {
        const std::scoped_lock lock{mutex_};
        return events_;
    }

    [[nodiscard]] auto cancellation_checks_after_request() const noexcept -> std::size_t {
        return cancellation_checks_after_request_.load(std::memory_order_relaxed);
    }

  private:
    mutable std::mutex mutex_;
    mutable std::vector<RecordedProgress> events_;
    mutable std::atomic<bool> cancel_{false};
    mutable std::atomic<std::size_t> cancellation_checks_after_request_{0};
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

class ThrowOnEveryCallbackObserver final : public calculation::CalculationObserver {
  public:
    void on_progress(const calculation::CalculationProgress& /*progress*/) const override {
        callbacks_.fetch_add(1, std::memory_order_relaxed);
        throw std::runtime_error{"observer failure"};
    }

    [[nodiscard]] auto callbacks() const noexcept -> std::size_t {
        return callbacks_.load(std::memory_order_relaxed);
    }

  private:
    mutable std::atomic<std::size_t> callbacks_{0};
};

class DirectTestMethod final : public methods::Method {
  public:
    explicit DirectTestMethod(const bool fails = false) : fails_{fails} {}

    [[nodiscard]] auto metadata() const noexcept -> const methods::MethodMetadata& override {
        return metadata_;
    }

    [[nodiscard]] auto requirements() const -> methods::MethodRequirements override {
        return {.coordinates = true,
                .resources = {.supports_cutoff = true,
                              .supports_cover = true,
                              .fragment_target_charge_policy =
                                  methods::FragmentTargetChargePolicy::zero}};
    }

    [[nodiscard]] auto option_schema() const noexcept
        -> std::span<const methods::MethodOptionSpec> override {
        return {};
    }

    [[nodiscard]] auto calculate(const methods::CalculationInput& input) const
        -> charges::AtomicCharges override {
        if (fails_) {
            throw std::logic_error{"direct observer test failure"};
        }
        return charges::AtomicCharges{std::vector<double>(input.molecule().atom_count())};
    }

  private:
    methods::MethodMetadata metadata_{.id = "direct-test",
                                      .name = "Direct test",
                                      .full_name = "Direct observer test method",
                                      .publication = std::nullopt,
                                      .priority = 0};
    bool fails_;
};

auto assert_single_terminal_fragment_progress(const std::vector<RecordedProgress>& events,
                                              const std::size_t fragment_count) -> void {
    REQUIRE(!events.empty());
    CHECK(std::count_if(events.begin(), events.end(), [fragment_count](const auto& event) {
              return event.completed_fragment_count == fragment_count &&
                     event.fragment_count == fragment_count;
          }) == 1);
    for (const auto& event : events) {
        CHECK(event.completed_fragment_count > 0);
        CHECK(event.completed_fragment_count <= event.fragment_count);
        CHECK(event.fragment_count == fragment_count);
    }
    for (std::size_t index = 1; index < events.size(); ++index) {
        CHECK(events[index - 1].completed_fragment_count < events[index].completed_fragment_count);
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

auto make_many_separated_waters() -> core::Molecule {
    auto atoms = std::vector<core::Atom>{};
    auto bonds = std::vector<core::Bond>{};
    auto positions = std::vector<core::Position>{};
    constexpr std::size_t water_count = 12;
    atoms.reserve(water_count * 3);
    bonds.reserve(water_count * 2);
    positions.reserve(water_count * 3);

    for (std::size_t water_index = 0; water_index < water_count; ++water_index) {
        const auto atom_index = water_index * 3;
        const auto x = static_cast<double>(water_index) * 10.0;
        atoms.insert(atoms.end(), {core::Atom{8}, core::Atom{1}, core::Atom{1}});
        bonds.emplace_back(atom_index, atom_index + 1, core::BondOrder::SINGLE);
        bonds.emplace_back(atom_index, atom_index + 2, core::BondOrder::SINGLE);
        positions.insert(positions.end(),
                         {{x, 0.0, 0.0}, {x + 0.96, 0.0, 0.0}, {x - 0.24, 0.93, 0.0}});
    }

    return core::Molecule{std::move(atoms),
                          std::move(bonds),
                          {core::Conformer{std::move(positions)}},
                          "many-separated-waters"};
}

auto assert_computation_boundary(const std::vector<RecordedProgress>& events,
                                 const calculation::ExecutionMode mode) -> void {
    REQUIRE(!events.empty());
    CHECK(events.front().phase == calculation::CalculationPhase::computation_started);
    CHECK(events.back().phase == calculation::CalculationPhase::computation_finished);
    CHECK(events.front().mode == mode);
    CHECK(events.back().mode == mode);
    CHECK(std::count_if(events.begin(), events.end(), [](const auto& event) {
              return event.phase == calculation::CalculationPhase::computation_started;
          }) == 1);
    CHECK(std::count_if(events.begin(), events.end(), [](const auto& event) {
              return event.phase == calculation::CalculationPhase::computation_finished;
          }) == 1);
}

} // namespace

TEST_CASE("default observer preserves successful calculation results", "[calculation][observer]") {
    {
        const auto result = calculate_application(calculation::AssessmentRequest{
            .molecules = core::MoleculeCollection{std::vector{chargefw::test::make_water()}},
            .parameter_sets = {},
            .method_id = "formal",
            .execution_selection =
                calculation::ExecutionSelection{calculation::ExecutionSelectionKind::full}});

        REQUIRE(result.calculated());
        CHECK(!result.cancelled());
        CHECK(result.charges->method_id() == std::string_view{"formal"});
    }
}

TEST_CASE("observer emits ordered computation and target phases", "[calculation][observer]") {
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

        REQUIRE(result.calculated());
        CHECK(!result.cancelled());

        const auto events = observer.events();
        REQUIRE(!events.empty());

        // The computation phase starts observation.
        CHECK(events[0].phase == calculation::CalculationPhase::computation_started);

        // computation_started must carry mode and method_id.
        CHECK(events[0].mode == calculation::ExecutionMode::full);
        CHECK(events[0].method_id == std::string_view{"formal"});

        // computation_finished must be the last event and carry mode/method_id.
        CHECK(events.back().phase == calculation::CalculationPhase::computation_finished);
        CHECK(events.back().mode == calculation::ExecutionMode::full);
        CHECK(events.back().method_id == std::string_view{"formal"});
        CHECK(std::count_if(events.begin(), events.end(), [](const auto& event) {
                  return event.phase == calculation::CalculationPhase::computation_started;
              }) == 1);
        CHECK(std::count_if(events.begin(), events.end(), [](const auto& event) {
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
        CHECK(target_started_it != computation_end_it);

        const auto target_finished_it =
            std::find_if(target_started_it, computation_end_it, [](const auto& e) {
                return e.phase == calculation::CalculationPhase::target_finished;
            });
        CHECK(target_finished_it != computation_end_it);
    }
}

TEST_CASE("assessment remains outside calculation observation", "[calculation][observer]") {
    {
        const auto observer = RecordingObserver{};
        const auto assessment = calculation::assess(calculation::AssessmentRequest{
            .molecules = core::MoleculeCollection{std::vector{chargefw::test::make_water()}},
            .parameter_sets = {},
            .method_id = "formal",
            .execution_selection =
                calculation::ExecutionSelection{calculation::ExecutionSelectionKind::full}});

        CHECK(assessment.executable());
        CHECK(observer.events().empty());
    }
}

TEST_CASE("resource threshold warnings remain assessment data before computation",
          "[calculation][observer]") {
    {
        const auto observer = RecordingObserver{};
        auto assessment = calculation::assess(calculation::AssessmentRequest{
            .molecules = core::MoleculeCollection{std::vector{chargefw::test::make_water()}},
            .parameter_sets = {},
            .method_id = "mgc",
            .execution_selection =
                calculation::ExecutionSelection{calculation::ExecutionSelectionKind::full},
            .resource_policy = {.cutoff_atom_threshold = 2}});

        REQUIRE(assessment.execution_issues().size() == 1);
        CHECK(assessment.execution_issues()[0].kind ==
              methods::ExecutionIssueKind::resource_threshold_exceeded);
        CHECK(observer.events().empty());

        const auto result = calculation::calculate(std::move(assessment), 1, observer);
        REQUIRE(result.calculated());
        CHECK(!result.cancelled());
        const auto events = observer.events();
        REQUIRE(!events.empty());
        CHECK(events[0].phase == calculation::CalculationPhase::computation_started);
    }
}

TEST_CASE("cancellation produces a terminal observer event", "[calculation][observer]") {
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

        CHECK(result.cancelled());
        CHECK(!result.calculated());
        CHECK(result.status == calculation::ExecutionStatus::cancelled);
        CHECK(result.selected_method_id == "formal");
        const auto events = observer.events();
        REQUIRE(!events.empty());
        CHECK(events.front().phase == calculation::CalculationPhase::computation_started);
        CHECK(events.back().phase == calculation::CalculationPhase::computation_finished);
        CHECK(std::count_if(events.begin(), events.end(), [](const auto& event) {
                  return event.phase == calculation::CalculationPhase::computation_finished;
              }) == 1);
    }
}

TEST_CASE("validation failures finish observation and propagate unchanged",
          "[calculation][observer]") {
    {
        const auto observer = RecordingObserver{};
        auto assessment = calculation::assess(calculation::AssessmentRequest{
            .molecules = core::MoleculeCollection{std::vector{chargefw::test::make_water()}},
            .parameter_sets = {},
            .method_id = "formal",
            .execution_selection =
                calculation::ExecutionSelection{calculation::ExecutionSelectionKind::full}});

        const auto calculate_with_invalid_thread_count = [&] -> void {
            static_cast<void>(calculation::calculate(
                std::move(assessment), std::numeric_limits<std::size_t>::max(), observer));
        };
        CHECK_THROWS_AS(calculate_with_invalid_thread_count(), std::invalid_argument);

        const auto events = observer.events();
        REQUIRE(!events.empty());
        CHECK(events.front().phase == calculation::CalculationPhase::computation_started);
        CHECK(events.back().phase == calculation::CalculationPhase::computation_finished);
        CHECK(events.back().mode == calculation::ExecutionMode::full);
        CHECK(events.back().method_id == std::string_view{"formal"});
        CHECK(std::count_if(events.begin(), events.end(), [](const auto& event) {
                  return event.phase == calculation::CalculationPhase::computation_finished;
              }) == 1);
    }
}

TEST_CASE("solver failures finish observation with a numerical result", "[calculation][observer]") {
    {
        const auto observer = RecordingObserver{};
        auto assessment = calculation::assess(calculation::AssessmentRequest{
            .molecules = core::MoleculeCollection{std::vector{chargefw::test::make_water()}},
            .parameter_sets = {make_invalid_qeq_parameters()},
            .method_id = "qeq",
            .parameter_set_id = "invalid-qeq",
            .execution_selection =
                calculation::ExecutionSelection{calculation::ExecutionSelectionKind::full}});

        const auto result = calculation::calculate(std::move(assessment), 1, observer);
        CHECK(result.status == calculation::ExecutionStatus::numerical_failure);
        CHECK_FALSE(result.calculated());
        REQUIRE(result.failure_message.has_value());
        CHECK(result.failure_message->contains("method 'qeq', molecule 1 ('water')"));

        const auto events = observer.events();
        REQUIRE(!events.empty());
        CHECK(events.front().phase == calculation::CalculationPhase::computation_started);
        CHECK(events.back().phase == calculation::CalculationPhase::computation_finished);
        CHECK(events.back().mode == calculation::ExecutionMode::full);
        CHECK(events.back().method_id == std::string_view{"qeq"});
        CHECK(std::count_if(events.begin(), events.end(), [](const auto& event) {
                  return event.phase == calculation::CalculationPhase::computation_finished;
              }) == 1);
    }
}

TEST_CASE("terminal observer callback failures do not terminate calculation",
          "[calculation][observer]") {
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

        REQUIRE(result.calculated());
        CHECK(!result.cancelled());
        CHECK(observer.terminal_callbacks() == 1);
    }
}

TEST_CASE("serial cutoff execution emits aggregate fragment progress", "[calculation][observer]") {
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

        REQUIRE(result.calculated());
        CHECK(!result.cancelled());

        const auto fragment_events = fragment_progress_events(observer);
        REQUIRE(!fragment_events.empty());
        assert_single_terminal_fragment_progress(fragment_events,
                                                 fragment_events.front().fragment_count);
    }
}

TEST_CASE("serial cover execution emits aggregate multi-pivot progress",
          "[calculation][observer]") {
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

        REQUIRE(result.calculated());
        CHECK(!result.cancelled());

        const auto fragment_events = fragment_progress_events(observer);
        REQUIRE(!fragment_events.empty());
        assert_single_terminal_fragment_progress(fragment_events,
                                                 fragment_events.front().fragment_count);
    }
}

TEST_CASE("multi-molecule target events carry source molecule identity",
          "[calculation][observer]") {
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

        REQUIRE(result.calculated());
        CHECK(!result.cancelled());

        const auto events = observer.events();

        // Collect all target_started events.
        std::vector<std::size_t> molecule_indices;
        for (const auto& event : events) {
            if (event.phase == calculation::CalculationPhase::target_started) {
                molecule_indices.push_back(event.molecule_index);
            }
        }

        REQUIRE(molecule_indices.size() == 2);
        CHECK(molecule_indices[0] == 0);
        CHECK(molecule_indices[1] == 1);

        // target_count must be consistent across events.
        for (const auto& event : events) {
            if (event.phase == calculation::CalculationPhase::target_started ||
                event.phase == calculation::CalculationPhase::target_finished) {
                CHECK(event.target_count == 2);
            }
        }
    }
}

TEST_CASE("empty reduced targets do not emit fragment progress", "[calculation][observer]") {
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
        REQUIRE(result.calculated());
        CHECK(fragment_progress_events(empty_observer).empty());
    }
}

TEST_CASE("parallel reduced targets each emit one terminal progress snapshot",
          "[calculation][observer]") {
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
        CHECK(result.calculated());

        for (const auto target_index : {std::size_t{0}, std::size_t{1}}) {
            auto target_events = std::vector<RecordedProgress>{};
            for (const auto& event : fragment_progress_events(observer)) {
                if (event.target_index == target_index) {
                    target_events.push_back(event);
                }
            }
            REQUIRE(!target_events.empty());
            assert_single_terminal_fragment_progress(target_events,
                                                     target_events.front().fragment_count);
        }
    }
}

TEST_CASE("observer target events preserve source target identity in every execution mode",
          "[calculation][observer]") {
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
        REQUIRE(result.calculated());
        REQUIRE(result.charges->size() == 3);

        const auto expected_targets = std::vector<charges::ChargeTarget>{
            {.molecule_index = 0, .conformer_index = 0},
            {.molecule_index = 0, .conformer_index = 1},
            {.molecule_index = 1, .conformer_index = 0},
        };
        for (std::size_t index = 0; index < expected_targets.size(); ++index) {
            CHECK(result.charges->assignment(index).target.molecule_index ==
                  expected_targets[index].molecule_index);
            CHECK(result.charges->assignment(index).target.conformer_index ==
                  expected_targets[index].conformer_index);
        }

        auto target_starts = std::vector<RecordedProgress>{};
        for (const auto& event : observer.events()) {
            if (event.phase == calculation::CalculationPhase::target_started) {
                target_starts.push_back(event);
            }
        }
        REQUIRE(target_starts.size() == expected_targets.size());
        for (std::size_t index = 0; index < target_starts.size(); ++index) {
            CHECK(target_starts[index].target_index == index);
            CHECK(target_starts[index].target_count == expected_targets.size());
            CHECK(target_starts[index].molecule_index == expected_targets[index].molecule_index);
            CHECK(target_starts[index].conformer_index == expected_targets[index].conformer_index);
        }
    }
}

TEST_CASE("no-plan failures occur before calculation observation begins",
          "[calculation][observer]") {
    {
        const auto observer = RecordingObserver{};
        auto assessment = calculation::assess(calculation::AssessmentRequest{
            .molecules = core::MoleculeCollection{std::vector{chargefw::test::make_water()}},
            .method_id = "formal",
            .execution_selection = calculation::ExecutionSelection{
                calculation::ExecutionSelectionKind::cover, calculation::minimum_reduced_radius}});
        CHECK(!assessment.executable());
        const auto result = calculation::calculate(std::move(assessment), 1, observer);
        CHECK(result.status == calculation::ExecutionStatus::no_executable_plan);
        CHECK_FALSE(result.calculated());
        CHECK(observer.events().empty());
    }
}

TEST_CASE("reduced fragment failures finish observation with a numerical result",
          "[calculation][observer]") {
    for (const auto mode : {calculation::ExecutionSelectionKind::cutoff,
                            calculation::ExecutionSelectionKind::cover}) {
        const auto observer = RecordingObserver{};
        auto assessment = calculation::assess(calculation::AssessmentRequest{
            .molecules = core::MoleculeCollection{std::vector{chargefw::test::make_water()}},
            .parameter_sets = {make_invalid_qeq_parameters()},
            .method_id = "qeq",
            .parameter_set_id = "invalid-qeq",
            .execution_selection =
                calculation::ExecutionSelection{mode, calculation::minimum_reduced_radius},
            .resource_policy = {.max_threads = 1}});

        const auto result = calculation::calculate(std::move(assessment), 1, observer);
        CHECK(result.status == calculation::ExecutionStatus::numerical_failure);
        CHECK_FALSE(result.calculated());
        REQUIRE(result.failure_message.has_value());
        CHECK(result.failure_message->contains("method 'qeq', molecule 1 ('water')"));

        const auto events = observer.events();
        REQUIRE(!events.empty());
        CHECK(events.front().phase == calculation::CalculationPhase::computation_started);
        CHECK(events.back().phase == calculation::CalculationPhase::computation_finished);
        CHECK(std::count_if(events.begin(), events.end(), [](const auto& event) {
                  return event.phase == calculation::CalculationPhase::computation_finished;
              }) == 1);
    }
}

TEST_CASE("observer callback failures do not alter calculation control flow",
          "[calculation][observer]") {
    {
        const auto observer = ThrowOnEveryCallbackObserver{};
        const auto result = calculate_application(
            calculation::AssessmentRequest{
                .molecules = core::MoleculeCollection{std::vector{chargefw::test::make_water()}},
                .parameter_sets = {},
                .method_id = "formal",
                .execution_selection =
                    calculation::ExecutionSelection{calculation::ExecutionSelectionKind::full}},
            observer);
        REQUIRE(result.calculated());
        CHECK(!result.cancelled());
        CHECK(observer.callbacks() >= 4);
    }
}

TEST_CASE("reduced execution observes cancellation after fragment progress",
          "[calculation][observer]") {
    for (const auto mode : {calculation::ExecutionSelectionKind::cutoff,
                            calculation::ExecutionSelectionKind::cover}) {
        for (const auto max_threads : {std::size_t{1}, std::size_t{2}}) {
            const auto observer = CancelAfterFirstFragmentProgress{};
            const auto result = calculate_application(
                calculation::AssessmentRequest{
                    .molecules =
                        core::MoleculeCollection{std::vector{make_many_separated_waters()}},
                    .parameter_sets = {},
                    .execution_selection =
                        calculation::ExecutionSelection{mode, calculation::minimum_reduced_radius},
                    .resource_policy = {.max_threads = max_threads}},
                observer);

            CHECK(result.cancelled());
            CHECK(!result.calculated());
            CHECK(result.status == calculation::ExecutionStatus::cancelled);
            CHECK(observer.cancellation_checks_after_request() > 0);

            const auto events = observer.events();
            REQUIRE(!events.empty());
            CHECK(std::any_of(events.begin(), events.end(), [](const auto& event) {
                return event.phase == calculation::CalculationPhase::fragment_progress;
            }));
            CHECK(std::none_of(events.begin(), events.end(), [](const auto& event) {
                return event.phase == calculation::CalculationPhase::target_finished;
            }));
            CHECK(std::count_if(events.begin(), events.end(), [](const auto& event) {
                      return event.phase == calculation::CalculationPhase::computation_finished;
                  }) == 1);
            CHECK(events.back().phase == calculation::CalculationPhase::computation_finished);
        }
    }
}

TEST_CASE("direct calculation emits terminal observation boundaries in every mode",
          "[calculation][observer]") {
    for (const auto mode : {calculation::ExecutionMode::full, calculation::ExecutionMode::cutoff,
                            calculation::ExecutionMode::cover}) {
        const auto observer = RecordingObserver{};
        const auto molecules = core::MoleculeCollection{std::vector{chargefw::test::make_water()}};
        const auto prepared = features::PreparedMoleculeCollection{molecules};
        const auto method = DirectTestMethod{};
        const auto selected = methods::ApplicableMethod{.method = &method};
        const auto policy =
            mode == calculation::ExecutionMode::full
                ? calculation::ExecutionPolicy{}
                : calculation::ExecutionPolicy{mode, calculation::minimum_reduced_radius,
                                               calculation::ChargeCorrectionPolicy::uniform};

        const auto result = calculation::calculate({.molecules = prepared,
                                                    .selected = selected,
                                                    .execution_policy = policy,
                                                    .max_threads = 1,
                                                    .observer = observer});
        CHECK(result.charges.size() == 1);
        assert_computation_boundary(observer.events(), mode);
    }
}

TEST_CASE("direct failures and cancellation finish observation in every mode",
          "[calculation][observer]") {
    for (const auto mode : {calculation::ExecutionMode::full, calculation::ExecutionMode::cutoff,
                            calculation::ExecutionMode::cover}) {
        const auto policy =
            mode == calculation::ExecutionMode::full
                ? calculation::ExecutionPolicy{}
                : calculation::ExecutionPolicy{mode, calculation::minimum_reduced_radius,
                                               calculation::ChargeCorrectionPolicy::uniform};

        {
            const auto observer = RecordingObserver{};
            const auto molecules =
                core::MoleculeCollection{std::vector{chargefw::test::make_water()}};
            const auto prepared = features::PreparedMoleculeCollection{molecules};
            const auto method = DirectTestMethod{true};
            const auto selected = methods::ApplicableMethod{.method = &method};
            const auto calculate_direct_failure = [&] -> void {
                static_cast<void>(calculation::calculate({.molecules = prepared,
                                                          .selected = selected,
                                                          .execution_policy = policy,
                                                          .max_threads = 1,
                                                          .observer = observer}));
            };
            CHECK_THROWS_AS(calculate_direct_failure(), std::exception);
            assert_computation_boundary(observer.events(), mode);
        }

        {
            const auto observer = CancelAfterFirstTarget{};
            const auto molecules =
                core::MoleculeCollection{std::vector{chargefw::test::make_water()}};
            const auto prepared = features::PreparedMoleculeCollection{molecules};
            const auto method = DirectTestMethod{};
            const auto selected = methods::ApplicableMethod{.method = &method};
            const auto calculate_direct_cancellation = [&] -> void {
                static_cast<void>(calculation::calculate({.molecules = prepared,
                                                          .selected = selected,
                                                          .execution_policy = policy,
                                                          .max_threads = 1,
                                                          .observer = observer}));
            };
            CHECK_THROWS_AS(calculate_direct_cancellation(), calculation::CalculationCancelled);
            assert_computation_boundary(observer.events(), mode);
        }
    }
}
