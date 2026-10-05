#include "calculation/fixed_charge_partition.h"
#include "calculation/full_execution.h"

#include "support/test_parameters.h"

#include <chargefw/core/atom.h>
#include <chargefw/core/bond.h>
#include <chargefw/core/conformer.h>
#include <chargefw/core/molecule.h>
#include <chargefw/core/molecule_collection.h>
#include <chargefw/features/prepared_molecule_collection.h>
#include <chargefw/methods/method.h>
#include <chargefw/methods/method_metadata.h>
#include <chargefw/methods/method_options.h>
#include <chargefw/methods/method_registry.h>
#include <chargefw/methods/method_requirements.h>
#include <chargefw/parameters/classification/parameter_classification.h>
#include <chargefw/parameters/models/atom_parameters.h>
#include <chargefw/parameters/models/bond_parameters.h>
#include <chargefw/parameters/models/common_parameters.h>
#include <chargefw/parameters/models/parameter_set.h>
#include <chargefw/parameters/models/parameter_set_metadata.h>

#include <snitch/snitch.hpp>

#include <array>
#include <cmath>
#include <optional>
#include <ranges>
#include <span>
#include <stdexcept>
#include <string>
#include <string_view>
#include <utility>
#include <vector>

namespace calculation = chargefw::calculation;
namespace charges = chargefw::charges;
namespace core = chargefw::core;
namespace features = chargefw::features;
namespace methods = chargefw::methods;

namespace {

class CapturingMethod final : public methods::Method {
  public:
    explicit CapturingMethod(const bool supports_fixed_point_sources)
        : supports_fixed_point_sources_{supports_fixed_point_sources} {}

    [[nodiscard]] auto metadata() const noexcept -> const methods::MethodMetadata& override {
        static constexpr methods::MethodMetadata value{.id = "capture",
                                                       .name = "Capture",
                                                       .full_name = "Capture",
                                                       .publication = std::nullopt,
                                                       .priority = 0};
        return value;
    }

    [[nodiscard]] auto requirements() const -> methods::MethodRequirements override {
        auto value = methods::MethodRequirements{};
        value.coordinates = true;
        value.supports_fixed_point_sources = supports_fixed_point_sources_;
        return value;
    }

    [[nodiscard]] auto option_schema() const noexcept
        -> std::span<const methods::MethodOptionSpec> override {
        return {};
    }

    [[nodiscard]] auto calculate(const methods::CalculationInput& input) const
        -> charges::AtomicCharges override {
        target_charges.push_back(input.target_charge());
        auto sources = std::vector<methods::FixedPointSource>{input.fixed_sources().begin(),
                                                              input.fixed_sources().end()};
        observed_sources.push_back(std::move(sources));
        return charges::AtomicCharges{
            std::vector<double>(input.molecule().atom_count(), input.target_charge())};
    }

    mutable std::vector<double> target_charges;
    mutable std::vector<std::vector<methods::FixedPointSource>> observed_sources;

  private:
    bool supports_fixed_point_sources_;
};

auto make_execution_partition() -> calculation::detail::FixedChargePartition {
    const auto molecules = core::MoleculeCollection{
        std::vector{core::Molecule{std::vector{core::Atom{6, 1}, core::Atom{12, 2}},
                                   {},
                                   {core::Conformer{{{0.0, 0.0, 0.0}, {2.0, 0.0, 0.0}}, "a"},
                                    core::Conformer{{{0.0, 0.0, 1.0}, {2.0, 0.0, 1.0}}, "b"}},
                                   "affected"},
                    core::Molecule{std::vector{core::Atom{6, -1}},
                                   {},
                                   {core::Conformer{{{5.0, 0.0, 0.0}}, "identity"}},
                                   "unaffected"}}};
    return calculation::detail::make_fixed_charge_partition(molecules,
                                                            calculation::FixedIons{{{0, 1, 0.25}}});
}

} // namespace

TEST_CASE("full execution supplies partition budgets and conformer-local sources",
          "[calculation][fixed-charge-execution]") {
    const auto partition = make_execution_partition();
    const features::PreparedMoleculeCollection prepared{partition.active_molecules};
    const CapturingMethod method{true};
    const methods::ApplicableMethod selected{.method = &method, .parameter_set = nullptr};

    const auto serial = calculation::calculate_full_charges(
        selected, prepared, 1, calculation::default_calculation_observer(), &partition);
    REQUIRE(serial.size() == 3);
    CHECK(method.target_charges == std::vector<double>{1.0, 1.0, -1.0});
    REQUIRE(method.observed_sources.size() == 3);
    REQUIRE(method.observed_sources[0].size() == 1);
    CHECK(method.observed_sources[0][0].charge == 0.25);
    CHECK(method.observed_sources[0][0].position.z == 0.0);
    CHECK(method.observed_sources[1][0].position.z == 1.0);
    CHECK(method.observed_sources[2].empty());
    CHECK(serial.assignment(0).charges.size() == 1);
    CHECK(serial.assignment(0).charges[0] == 1.0);

    CHECK(partition.targets[0].sources[0].charge == 0.25);
}

TEST_CASE("full execution rejects unsupported fixed ions and mismatched prepared ownership",
          "[calculation][fixed-charge-execution]") {
    const auto partition = make_execution_partition();
    const features::PreparedMoleculeCollection prepared{partition.active_molecules};
    const CapturingMethod unsupported{false};
    const methods::ApplicableMethod unsupported_selected{.method = &unsupported,
                                                         .parameter_set = nullptr};
    CHECK_THROWS_AS(calculation::calculate_full_charges(unsupported_selected, prepared, 1,
                                                        calculation::default_calculation_observer(),
                                                        &partition),
                    std::invalid_argument);

    const core::MoleculeCollection other_molecules{partition.active_molecules};
    const features::PreparedMoleculeCollection mismatched{other_molecules};
    const CapturingMethod supported{true};
    const methods::ApplicableMethod selected{.method = &supported, .parameter_set = nullptr};
    CHECK_THROWS_AS(calculation::calculate_full_charges(selected, mismatched, 1,
                                                        calculation::default_calculation_observer(),
                                                        &partition),
                    std::invalid_argument);
}

TEST_CASE("full execution keeps ordinary no-partition behavior",
          "[calculation][fixed-charge-execution]") {
    const auto molecules = core::MoleculeCollection{std::vector{core::Molecule{
        std::vector{core::Atom{6, 3}}, {}, {core::Conformer{{{0.0, 0.0, 0.0}}}}, "ordinary"}}};
    const features::PreparedMoleculeCollection prepared{molecules};
    const CapturingMethod method{false};
    const methods::ApplicableMethod selected{.method = &method, .parameter_set = nullptr};
    const auto result = calculation::calculate_full_charges(
        selected, prepared, 1, calculation::default_calculation_observer());
    REQUIRE(result.size() == 1);
    CHECK(method.target_charges == std::vector<double>{3.0});
    REQUIRE(method.observed_sources.size() == 1);
    CHECK(method.observed_sources[0].empty());
}

TEST_CASE("parameterized full EEM uses partition sources and active budgets",
          "[calculation][fixed-charge-execution]") {
    constexpr auto kappa = 2.5;
    const auto original = core::MoleculeCollection{std::vector{
        core::Molecule{
            std::vector{core::Atom{1, 1}, core::Atom{8, -1}, core::Atom{12, 2}},
            {},
            {core::Conformer{{{0.0, 0.0, 0.0}, {2.0, 0.0, 0.0}, {0.0, 3.0, 0.0}}, "first"},
             core::Conformer{{{0.0, 0.0, 0.0}, {2.0, 0.0, 0.0}, {4.0, 1.0, 0.0}}, "second"}},
            "active-pair-and-source"},
        core::Molecule{std::vector{core::Atom{1, -1}},
                       {},
                       {core::Conformer{{{8.0, 0.0, 0.0}}, "unaffected"}},
                       "unaffected"}}};
    const auto partition = calculation::detail::make_fixed_charge_partition(
        original, calculation::FixedIons{{{0, 2, 0.4}}});
    const features::PreparedMoleculeCollection prepared{partition.active_molecules};
    const auto parameter_set = chargefw::parameters::ParameterSet{
        chargefw::parameters::ParameterSetMetadata{
            .id = "execution-eem", .method_id = "eem", .name = "Execution EEM"},
        chargefw::parameters::CommonParameters{{{.name = "kappa", .value = kappa}}},
        chargefw::parameters::AtomParameters{
            {{.key = chargefw::test::plain_atom_key(1),
              .parameters = {{.name = "A", .value = 1.0}, {.name = "B", .value = 5.0}}},
             {.key = chargefw::test::plain_atom_key(8),
              .parameters = {{.name = "A", .value = 2.0}, {.name = "B", .value = 9.0}}}}}};
    const auto* eem = methods::method_registry().find("eem");
    REQUIRE(eem != nullptr);
    const auto selected = methods::ApplicableMethod{
        .method = eem,
        .parameter_set = &parameter_set,
        .classifications = {chargefw::parameters::ParameterClassification{
                                chargefw::parameters::AtomParameterClassification{{0, 1}}},
                            chargefw::parameters::ParameterClassification{
                                chargefw::parameters::AtomParameterClassification{{0}}}}};

    const auto expected_pair = [](const double source_potential_first,
                                  const double source_potential_second) {
        constexpr auto active_charge = 0.0;
        constexpr auto cross_interaction = kappa / 2.0;
        const auto rhs_difference = -1.0 - source_potential_first + 2.0 + source_potential_second;
        const auto first = (rhs_difference + (9.0 - cross_interaction) * active_charge) /
                           (5.0 + 9.0 - 2.0 * cross_interaction);
        return std::vector<double>{first, active_charge - first};
    };
    const auto first_expected = expected_pair(kappa * 0.4 / 3.0, kappa * 0.4 / std::sqrt(13.0));
    const auto second_expected =
        expected_pair(kappa * 0.4 / std::sqrt(17.0), kappa * 0.4 / std::sqrt(5.0));

    const auto serial = calculation::calculate_full_charges(
        selected, prepared, 1, calculation::default_calculation_observer(), &partition);
    const auto parallel = calculation::calculate_full_charges(
        selected, prepared, 2, calculation::default_calculation_observer(), &partition);
    const auto repeated = calculation::calculate_full_charges(
        selected, prepared, 1, calculation::default_calculation_observer(), &partition);
    REQUIRE(serial.size() == 3);
    REQUIRE(parallel.size() == serial.size());
    REQUIRE(repeated.size() == serial.size());
    for (std::size_t index = 0; index < serial.size(); ++index) {
        CHECK(std::ranges::equal(serial.assignment(index).charges.values(),
                                 repeated.assignment(index).charges.values()));
        CHECK(std::ranges::equal(serial.assignment(index).charges.values(),
                                 parallel.assignment(index).charges.values()));
    }
    for (std::size_t conformer = 0; conformer < 2; ++conformer) {
        const auto& expected = conformer == 0 ? first_expected : second_expected;
        const auto& values = serial.assignment(conformer).charges;
        CHECK(std::abs(values[0] - expected[0]) < 1e-10);
        CHECK(std::abs(values[1] - expected[1]) < 1e-10);
        CHECK(std::abs(values.total()) < 1e-12);
    }
    CHECK(std::abs(serial.assignment(2).charges[0] + 1.0) < 1e-12);
    CHECK(partition.targets[0].sources[0].charge == 0.4);
    CHECK(partition.targets[0].source_positions[0][0].y == 3.0);
    CHECK(partition.targets[0].source_positions[1][0].x == 4.0);
    CHECK(partition.targets[1].sources.empty());
}

TEST_CASE("SQE family full-execution fixed Mg response decays with distance",
          "[calculation][fixed-charge-execution][sqe][sqeq0][sqeqp]") {
    constexpr auto distances = std::array{3.0, 12.0, 120.0};
    constexpr auto mg_charge = 2.0;
    auto conformers = std::vector<core::Conformer>{};
    for (const auto distance : distances) {
        conformers.emplace_back(
            std::vector<core::Position>{{0.0, 0.0, 0.0}, {1.5, 0.0, 0.0}, {-distance, 0.0, 0.0}});
    }
    const auto original = core::MoleculeCollection{std::vector{
        core::Molecule{std::vector{core::Atom{1, 1}, core::Atom{8, -1}, core::Atom{12, 2}},
                       {core::Bond{0, 1}},
                       std::move(conformers),
                       "active-pair-and-Mg"}}};
    const auto partition = calculation::detail::make_fixed_charge_partition(
        original, calculation::FixedIons{{{0, 2, mg_charge}}});
    const features::PreparedMoleculeCollection prepared{partition.active_molecules};

    for (const auto method_id : {"sqe", "sqeq0", "sqeqp"}) {
        CAPTURE(method_id);
        const auto parameter_set = chargefw::parameters::ParameterSet{
            chargefw::parameters::ParameterSetMetadata{.id =
                                                           std::string{"Mg-response-"} + method_id,
                                                       .method_id = method_id,
                                                       .name = "Mg response SQE parameters"},
            {},
            chargefw::parameters::AtomParameters{
                {{.key = chargefw::test::plain_atom_key(1),
                  .parameters = {{.name = "electronegativity", .value = 4.5280},
                                 {.name = "hardness", .value = 13.8904},
                                 {.name = "width", .value = 1.0},
                                 {.name = "q0", .value = 0.25}}},
                 {.key = chargefw::test::plain_atom_key(8),
                  .parameters = {{.name = "electronegativity", .value = 8.741},
                                 {.name = "hardness", .value = 13.364},
                                 {.name = "width", .value = 1.0},
                                 {.name = "q0", .value = -0.5}}}}},
            chargefw::parameters::BondParameters{
                {{.key = chargefw::test::single_bond_key(1, 8),
                  .parameters = {{.name = "kappa", .value = 1.0}}}}}};
        const auto* method = methods::method_registry().find(method_id);
        REQUIRE(method != nullptr);
        const auto selected = methods::ApplicableMethod{
            .method = method,
            .parameter_set = &parameter_set,
            .method_options = methods::make_default_options(method->option_schema()),
            .classifications = {chargefw::parameters::ParameterClassification{
                chargefw::parameters::AtomParameterClassification{{0, 1}},
                chargefw::parameters::BondParameterClassification{{0}}}}};
        const auto& observer = calculation::default_calculation_observer();
        // The ion-free reference reuses exactly the same prepared active geometry and charge.
        const auto ion_free = calculation::calculate_full_charges(selected, prepared, 1, observer);
        const auto full =
            calculation::calculate_full_charges(selected, prepared, 1, observer, &partition);
        REQUIRE(ion_free.size() == distances.size());
        REQUIRE(full.size() == distances.size());
        auto previous_response = 1.0;
        auto near_response = 0.0;
        for (std::size_t conformer = 0; conformer < distances.size(); ++conformer) {
            CAPTURE(distances[conformer]);
            const auto& reference = ion_free.assignment(conformer).charges;
            const auto& active = full.assignment(conformer).charges;
            REQUIRE(active.size() == 2);
            REQUIRE(reference.size() == active.size());
            CHECK(std::abs(reference.total() - partition.targets[0].active_charge) < 1e-12);
            CHECK(std::abs(active.total() - partition.targets[0].active_charge) < 1e-12);
            auto squared_response = 0.0;
            for (std::size_t atom = 0; atom < active.size(); ++atom) {
                const auto difference = active[atom] - reference[atom];
                squared_response += difference * difference;
            }
            const auto response = std::sqrt(squared_response);
            CHECK(response > 1e-8);
            if (conformer == 0) {
                near_response = response;
                CHECK(near_response > 1e-4);
            } else {
                CHECK(response < previous_response);
            }
            previous_response = response;
        }
        CHECK(previous_response < near_response * 0.01);
        CHECK(previous_response < 2e-5);
    }
}
