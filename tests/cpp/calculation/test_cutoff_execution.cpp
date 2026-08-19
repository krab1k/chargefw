#include "support/test_molecules.h"
#include "support/test_parameters.h"

#include <chargefw/calculation/calculation.h>
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
#include <chargefw/parameters/models/atom_parameters.h>
#include <chargefw/parameters/models/common_parameters.h>
#include <chargefw/parameters/models/parameter_key.h>
#include <chargefw/parameters/models/parameter_set.h>
#include <chargefw/parameters/models/parameter_set_metadata.h>

#include <cassert>
#include <cmath>
#include <optional>
#include <span>
#include <string_view>
#include <utility>
#include <vector>

namespace calculation = chargefw::calculation;
namespace charges = chargefw::charges;
namespace core = chargefw::core;
namespace features = chargefw::features;
namespace methods = chargefw::methods;
namespace parameters = chargefw::parameters;

namespace {

class ZeroFragmentMethod final : public methods::Method {
  public:
    [[nodiscard]] auto metadata() const noexcept -> const methods::MethodMetadata& override {
        static constexpr methods::MethodMetadata metadata{.id = "zero-fragment",
                                                          .name = "Zero fragment",
                                                          .full_name = "Zero fragment",
                                                          .publication = std::nullopt,
                                                          .priority = 0};
        return metadata;
    }

    [[nodiscard]] auto requirements() const -> methods::MethodRequirements override {
        auto requirements = methods::MethodRequirements{};
        requirements.coordinates = true;
        requirements.resources.supports_cutoff = true;
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
        return charges::AtomicCharges{std::vector<double>(input.molecule().atom_count(), 0.0)};
    }
};

auto make_eem_parameters() -> parameters::ParameterSet {
    return parameters::ParameterSet{
        parameters::ParameterSetMetadata{.id = "test-eem", .method_id = "eem", .name = "Test EEM"},
        parameters::CommonParameters{{{.name = "kappa", .value = 1.0}}},
        parameters::AtomParameters{
            {{.key = chargefw::test::plain_atom_key(1),
              .parameters = {{.name = "A", .value = 1.0}, {.name = "B", .value = 10.0}}},
             {.key = chargefw::test::plain_atom_key(8),
              .parameters = {{.name = "A", .value = 2.0}, {.name = "B", .value = 10.0}}}}}};
}

auto make_qeq_parameters() -> parameters::ParameterSet {
    return parameters::ParameterSet{
        parameters::ParameterSetMetadata{.id = "test-qeq", .method_id = "qeq", .name = "Test QEq"},
        {},
        parameters::AtomParameters{{{.key = chargefw::test::plain_atom_key(1),
                                     .parameters = {{.name = "electronegativity", .value = 4.5280},
                                                    {.name = "hardness", .value = 13.8904}}},
                                    {.key = chargefw::test::plain_atom_key(8),
                                     .parameters = {{.name = "electronegativity", .value = 8.741},
                                                    {.name = "hardness", .value = 13.364}}}}}};
}

auto assert_cutoff_matches_full(const std::string_view method_id,
                                std::vector<parameters::ParameterSet> parameter_sets = {}) -> void {
    const auto molecules =
        core::MoleculeCollection{std::vector{chargefw::test::make_two_conformer_water()}};
    const auto full = calculation::calculate(
        calculation::ApplicationCalculationRequest{.molecules = molecules,
                                                   .parameter_sets = parameter_sets,
                                                   .method_id = std::string{method_id}});
    const auto cutoff = calculation::calculate(calculation::ApplicationCalculationRequest{
        .molecules = molecules,
        .parameter_sets = std::move(parameter_sets),
        .method_id = std::string{method_id},
        .execution_selection =
            calculation::ExecutionSelection{calculation::ExecutionSelectionKind::cutoff, 8.0}});

    assert(full.calculated());
    assert(cutoff.calculated());
    assert(cutoff.execution_policy.has_value());
    assert(cutoff.execution_policy->mode() == calculation::ExecutionMode::cutoff);
    assert(cutoff.execution_policy->charge_correction() ==
           calculation::ChargeCorrectionPolicy::uniform);
    assert(full.charges->size() == cutoff.charges->size());

    for (std::size_t assignment_index = 0; assignment_index < full.charges->size();
         ++assignment_index) {
        const auto& full_charges = full.charges->assignment(assignment_index).charges;
        const auto& cutoff_charges = cutoff.charges->assignment(assignment_index).charges;
        assert(full_charges.size() == cutoff_charges.size());
        for (std::size_t atom_index = 0; atom_index < full_charges.size(); ++atom_index) {
            assert(std::abs(full_charges[atom_index] - cutoff_charges[atom_index]) < 1.0e-10);
        }
    }
}

} // namespace

auto main() -> int {
    const ZeroFragmentMethod zero_method;
    const auto charged_molecule = core::Molecule{
        std::vector{core::Atom{1, 1}, core::Atom{1, 0}},
        {},
        std::vector{core::Conformer{{core::Position{.x = 0.0}, core::Position{.x = 20.0}}}},
        "charged"};
    const auto collection = core::MoleculeCollection{std::vector{charged_molecule}};
    const features::PreparedMoleculeCollection prepared{collection};
    const methods::ApplicableMethod selected{.method = &zero_method, .parameter_set = nullptr};

    const auto corrected =
        calculation::calculate({.molecules = prepared,
                                .selected = selected,
                                .execution_policy = calculation::ExecutionPolicy{
                                    calculation::ExecutionMode::cutoff, 8.0,
                                    calculation::ChargeCorrectionPolicy::uniform}});
    assert(corrected.charges.assignment(0).charges[0] == 0.5);
    assert(corrected.charges.assignment(0).charges[1] == 0.5);

    const auto uncorrected = calculation::calculate(
        {.molecules = prepared,
         .selected = selected,
         .execution_policy = calculation::ExecutionPolicy{
             calculation::ExecutionMode::cutoff, 8.0, calculation::ChargeCorrectionPolicy::none}});
    assert(uncorrected.charges.assignment(0).charges[0] == 0.0);
    assert(uncorrected.charges.assignment(0).charges[1] == 0.0);

    assert_cutoff_matches_full("eem", {make_eem_parameters()});
    assert_cutoff_matches_full("qeq", {make_qeq_parameters()});
    assert_cutoff_matches_full("eqeq");
    return 0;
}
