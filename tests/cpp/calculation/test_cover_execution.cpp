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

#include <cassert>
#include <optional>
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
        ++calls;
        return charges::AtomicCharges{std::vector<double>(
            input.molecule().atom_count(), static_cast<double>(input.molecule().atom_count()))};
    }

    mutable std::size_t calls = 0;
};

auto make_linear_molecule() -> core::Molecule {
    auto atoms = std::vector<core::Atom>{};
    auto positions = std::vector<core::Position>{};
    for (std::size_t atom_index = 0; atom_index < 8; ++atom_index) {
        atoms.emplace_back(6);
        positions.push_back(core::Position{.x = static_cast<double>(atom_index)});
    }
    return core::Molecule{
        std::move(atoms), {}, std::vector{core::Conformer{std::move(positions)}}, "cover-line"};
}

} // namespace

auto main() -> int {
    const FragmentSizeMethod method;
    const auto collection = core::MoleculeCollection{std::vector{make_linear_molecule()}};
    const features::PreparedMoleculeCollection prepared{collection};
    const methods::ApplicableMethod selected{.method = &method, .parameter_set = nullptr};

    const auto result = calculation::calculate_cover_charges(
        selected, prepared,
        calculation::ExecutionPolicy{calculation::ExecutionMode::cover, 8.0,
                                     calculation::ChargeCorrectionPolicy::none});

    assert(method.calls == 2);
    assert(result.size() == 1);
    const auto& values = result.assignment(0).charges;
    assert(values.size() == 8);
    for (const auto value : values.values()) {
        assert(value == 8.0);
    }

    return 0;
}
