#include "calculation/cover_execution.h"

#include <chargefw/calculation/execution_policy.h>
#include <chargefw/charges/atomic_charges.h>
#include <chargefw/core/atom.h>
#include <chargefw/core/conformer.h>
#include <chargefw/core/molecule.h>
#include <chargefw/core/molecule_collection.h>
#include <chargefw/core/position.h>
#include <chargefw/features/prepared_molecule_collection.h>
#include <chargefw/methods/method.h>
#include <chargefw/methods/method_metadata.h>
#include <chargefw/methods/method_options.h>
#include <chargefw/methods/method_requirements.h>

#include <atomic>
#include <cmath>
#include <optional>
#include <snitch/snitch.hpp>
#include <span>
#include <string_view>
#include <vector>

namespace calculation = chargefw::calculation;
namespace charges = chargefw::charges;
namespace core = chargefw::core;
namespace features = chargefw::features;
namespace methods = chargefw::methods;

namespace {

class FragmentSizeMethod final : public methods::Method {
  public:
    [[nodiscard]] auto metadata() const noexcept -> const methods::MethodMetadata& override {
        static constexpr methods::MethodMetadata metadata{.id = "cover-fragment-size",
                                                          .name = "Cover fragment size",
                                                          .full_name = "Cover fragment size",
                                                          .publication = std::nullopt,
                                                          .priority = 0};
        return metadata;
    }

    [[nodiscard]] auto requirements() const -> methods::MethodRequirements override {
        auto requirements = methods::MethodRequirements{};
        requirements.coordinates = true;
        requirements.resources.supports_cover = true;
        requirements.resources.fragment_target_charge_policy =
            methods::FragmentTargetChargePolicy::proportional_to_atom_count;
        return requirements;
    }

    [[nodiscard]] auto option_schema() const noexcept
        -> std::span<const methods::MethodOptionSpec> override {
        return {};
    }

    [[nodiscard]] auto calculate(const methods::CalculationInput& input) const
        -> charges::AtomicCharges override {
        calls.fetch_add(1, std::memory_order_relaxed);
        return charges::AtomicCharges{std::vector<double>(
            input.molecule().atom_count(), static_cast<double>(input.molecule().atom_count()))};
    }

    mutable std::atomic_size_t calls = 0;
};

class TargetChargeMethod final : public methods::Method {
  public:
    [[nodiscard]] auto metadata() const noexcept -> const methods::MethodMetadata& override {
        static constexpr methods::MethodMetadata metadata{.id = "cover-target-charge",
                                                          .name = "Cover target charge",
                                                          .full_name = "Cover target charge",
                                                          .publication = std::nullopt,
                                                          .priority = 0};
        return metadata;
    }

    [[nodiscard]] auto requirements() const -> methods::MethodRequirements override {
        auto requirements = methods::MethodRequirements{};
        requirements.coordinates = true;
        requirements.resources.supports_cover = true;
        requirements.resources.fragment_target_charge_policy =
            methods::FragmentTargetChargePolicy::proportional_to_atom_count;
        return requirements;
    }

    [[nodiscard]] auto option_schema() const noexcept
        -> std::span<const methods::MethodOptionSpec> override {
        return {};
    }

    [[nodiscard]] auto calculate(const methods::CalculationInput& input) const
        -> charges::AtomicCharges override {
        return charges::AtomicCharges{
            std::vector<double>(input.molecule().atom_count(), input.target_charge())};
    }
};

auto make_linear_molecule(const std::size_t atom_count = 8, const int formal_charge = 0)
    -> core::Molecule {
    auto atoms = std::vector<core::Atom>{};
    auto positions = std::vector<core::Position>{};
    for (std::size_t atom_index = 0; atom_index < atom_count; ++atom_index) {
        atoms.emplace_back(6, formal_charge);
        positions.push_back(core::Position{.x = static_cast<double>(atom_index)});
    }
    return core::Molecule{
        std::move(atoms), {}, std::vector{core::Conformer{std::move(positions)}}, "cover-line"};
}

} // namespace

TEST_CASE("cover execution produces fragment-size-dependent charges",
          "[calculation][cover-execution]") {
    const auto collection = core::MoleculeCollection{std::vector{make_linear_molecule(10)}};
    const features::PreparedMoleculeCollection prepared{collection};
    const auto policy = calculation::ExecutionPolicy{calculation::ExecutionMode::cover, 8.0,
                                                     calculation::ChargeCorrectionPolicy::none};

    {
        const FragmentSizeMethod method;
        const methods::ApplicableMethod selected{.method = &method, .parameter_set = nullptr};
        const auto serial = calculation::calculate_cover_charges(
            selected, prepared, policy, 1, calculation::default_calculation_observer());
        const auto parallel = calculation::calculate_cover_charges(
            selected, prepared, policy, 0, calculation::default_calculation_observer());

        CHECK(method.calls.load() == 6);
        CHECK(serial.size() == 1);
        CHECK(parallel.size() == 1);
        const auto& serial_values = serial.assignment(0).charges;
        const auto& parallel_values = parallel.assignment(0).charges;
        CHECK(serial_values.size() == 10);
        CHECK(parallel_values.size() == serial_values.size());
        for (std::size_t atom_index = 0; atom_index < serial_values.size(); ++atom_index) {
            CHECK(parallel_values[atom_index] == serial_values[atom_index]);
            const auto expected = atom_index < 4 ? 9.0 : 10.0;
            CHECK(serial_values[atom_index] == expected);
        }
    }

    {
        const auto charged_collection =
            core::MoleculeCollection{std::vector{make_linear_molecule(10, 1)}};
        const features::PreparedMoleculeCollection charged_prepared{charged_collection};
        const TargetChargeMethod method;
        const methods::ApplicableMethod selected{.method = &method, .parameter_set = nullptr};
        const auto corrected = calculation::calculate_cover_charges(
            selected, charged_prepared,
            calculation::ExecutionPolicy{calculation::ExecutionMode::cover, 8.0,
                                         calculation::ChargeCorrectionPolicy::uniform},
            1, calculation::default_calculation_observer());

        const auto& values = corrected.assignment(0).charges;
        for (std::size_t atom_index = 0; atom_index < values.size(); ++atom_index) {
            const auto expected = atom_index < 4 ? 0.4 : 1.4;
            CHECK(std::abs(values[atom_index] - expected) < 1.0e-12);
        }
    }
}
