#include "calculation/target_execution.h"

#include <chargefw/core/atom.h>
#include <chargefw/core/bond.h>
#include <chargefw/core/conformer.h>
#include <chargefw/core/molecule.h>
#include <chargefw/core/molecule_collection.h>
#include <chargefw/features/prepared_molecule_collection.h>
#include <chargefw/methods/method.h>
#include <chargefw/methods/method_metadata.h>
#include <chargefw/methods/method_options.h>
#include <chargefw/methods/method_requirements.h>

#include <cstddef>
#include <optional>
#include <snitch/snitch.hpp>
#include <span>
#include <vector>

namespace calculation = chargefw::calculation;
namespace charges = chargefw::charges;
namespace core = chargefw::core;
namespace features = chargefw::features;
namespace methods = chargefw::methods;

namespace {

class TestMethod final : public methods::Method {
  public:
    [[nodiscard]] auto metadata() const noexcept -> const methods::MethodMetadata& override {
        static constexpr methods::MethodMetadata metadata{.id = "target-execution",
                                                          .name = "Target execution",
                                                          .full_name = "Target execution",
                                                          .publication = std::nullopt,
                                                          .priority = 0};
        return metadata;
    }

    [[nodiscard]] auto requirements() const -> methods::MethodRequirements override {
        return {};
    }

    [[nodiscard]] auto option_schema() const noexcept
        -> std::span<const methods::MethodOptionSpec> override {
        return {};
    }

    [[nodiscard]] auto calculate(const methods::CalculationInput&) const
        -> charges::AtomicCharges override {
        return charges::AtomicCharges{std::vector<double>{}};
    }
};

auto make_molecules(const std::size_t count) -> core::MoleculeCollection {
    auto molecules = std::vector<core::Molecule>{};
    molecules.reserve(count);
    for (std::size_t index = 0; index < count; ++index) {
        molecules.emplace_back(std::vector{core::Atom{6}}, std::vector<core::Bond>{},
                               std::vector<core::Conformer>{}, "target");
    }
    return core::MoleculeCollection{std::move(molecules)};
}

auto assert_worker_budget(const calculation::detail::ParallelizationLevel parallelization,
                          const std::size_t requested_threads, const std::size_t target_count,
                          const std::size_t expected_fragment_threads) -> void {
    const auto molecules = make_molecules(target_count);
    const features::PreparedMoleculeCollection prepared{molecules};
    const TestMethod method;
    const methods::ApplicableMethod selected{.method = &method, .parameter_set = nullptr};
    auto received_fragment_threads = std::vector<std::size_t>(target_count);
    const auto mode = parallelization == calculation::detail::ParallelizationLevel::targets
                          ? calculation::ExecutionMode::full
                          : calculation::ExecutionMode::cutoff;

    const auto result = calculation::detail::execute_calculation_targets(
        selected, prepared, mode, false, parallelization, requested_threads,
        calculation::default_calculation_observer(),
        [&](const features::PreparedMolecule& molecule,
            const chargefw::parameters::ParameterClassification*, const std::optional<std::size_t>,
            const std::size_t fragment_threads,
            const calculation::detail::ProgressContext& context) {
            received_fragment_threads[context.target_index] = fragment_threads;
            return charges::AtomicCharges{
                std::vector<double>(molecule.molecule().atom_count(), 0.0)};
        });

    REQUIRE(result.size() == target_count);
    for (std::size_t target_index = 0; target_index < target_count; ++target_index) {
        CHECK(result.assignment(target_index).target.molecule_index == target_index);
        CHECK(received_fragment_threads[target_index] == expected_fragment_threads);
    }
}

} // namespace

TEST_CASE("calculation targets distribute worker budget across parallelization levels",
          "[calculation][calculation-targets]") {
    for (const auto requested_threads : {std::size_t{0}, std::size_t{1}, std::size_t{2}}) {
        for (const auto target_count : {std::size_t{0}, std::size_t{1}, std::size_t{2}}) {
            assert_worker_budget(calculation::detail::ParallelizationLevel::targets,
                                 requested_threads, target_count, 1);
            assert_worker_budget(calculation::detail::ParallelizationLevel::fragments,
                                 requested_threads, target_count, requested_threads);
        }
    }
}
