#include "support/test_calculation.h"
#include "support/test_molecules.h"
#include "support/test_parameters.h"

#include <chargefw/core/atom.h>
#include <chargefw/core/conformer.h>
#include <chargefw/core/molecule.h>
#include <chargefw/core/position.h>
#include <chargefw/features/conformer_features.h>
#include <chargefw/features/prepared_molecule.h>
#include <chargefw/methods/calculation_input.h>
#include <chargefw/methods/method_registry.h>
#include <chargefw/parameters/classification/parameter_classification.h>
#include <chargefw/parameters/models/atom_parameters.h>
#include <chargefw/parameters/models/bond_parameters.h>
#include <chargefw/parameters/models/parameter_key.h>
#include <chargefw/parameters/models/parameter_set.h>
#include <chargefw/parameters/models/parameter_set_metadata.h>
#include <chargefw/parameters/models/parameter_view.h>

#include <array>
#include <cmath>
#include <numbers>
#include <snitch/snitch.hpp>
#include <stdexcept>
#include <string>
#include <string_view>
#include <utility>
#include <vector>

namespace parameters = chargefw::parameters;

namespace {

auto calculate_pair(const chargefw::core::Molecule& molecule, const std::string_view method_id,
                    const double target_charge,
                    const std::span<const chargefw::methods::FixedPointSource> sources)
    -> chargefw::charges::AtomicCharges {
    const chargefw::features::PreparedMolecule prepared{molecule};
    const chargefw::features::ConformerFeatures geometry{molecule, 0};
    const auto parameter_set = chargefw::test::make_sqe_ho_parameters(method_id);
    const chargefw::parameters::ParameterClassification classification{
        chargefw::parameters::AtomParameterClassification{std::vector<std::size_t>{0, 1}},
        chargefw::parameters::BondParameterClassification{std::vector<std::size_t>{0}}};
    const chargefw::parameters::ParameterView parameters{parameter_set, classification};
    const auto* method = chargefw::methods::method_registry().find(method_id);
    if (method == nullptr) {
        throw std::logic_error{"SQE method is not registered"};
    }
    const auto options = chargefw::methods::make_default_options(method->option_schema());
    const chargefw::methods::CalculationInput input{prepared,  options,     target_charge,
                                                    &geometry, &parameters, sources};
    return method->calculate(input);
}

} // namespace

TEST_CASE("SQE fixed ion right-hand sides preserve conserved pair totals",
          "[methods][sqe][sqeq0]") {
    using namespace chargefw;
    const auto make_pair = [](const int hydrogen_charge, const int oxygen_charge) {
        return core::Molecule{
            std::vector{core::Atom{1, hydrogen_charge}, core::Atom{8, oxygen_charge}},
            {core::Bond{0, 1}},
            {core::Conformer{{core::Position{}, core::Position{.x = 1.0}}}},
            "fixed-ion-pair"};
    };
    const std::array source{
        methods::FixedPointSource{.position = core::Position{.x = 2.0}, .charge = 0.5}};
    constexpr auto hydrogen_electronegativity = 4.528;
    constexpr auto oxygen_electronegativity = 8.741;
    constexpr auto hydrogen_hardness = 13.8904;
    constexpr auto oxygen_hardness = 13.364;
    constexpr auto kappa = 1.0;
    const auto bond_interaction = std::erf(1.0 / 2.0) / 1.0;
    const auto source_interaction = [](const double distance) {
        return std::erf(distance / std::numbers::sqrt2) / distance;
    };
    const auto source_h = 0.5 * source_interaction(2.0);
    const auto source_o = 0.5 * source_interaction(1.0);
    const auto denominator = hydrogen_hardness + oxygen_hardness - 2.0 * bond_interaction + kappa;

    const auto sqe_transfer =
        ((-hydrogen_electronegativity - source_h) - (-oxygen_electronegativity - source_o)) /
        denominator;
    const auto sqe = calculate_pair(make_pair(0, 0), "sqe", 0.0, source);
    CHECK(std::abs(sqe.values()[0] - sqe_transfer) < 1.0e-12);
    CHECK(std::abs(sqe.values()[1] + sqe_transfer) < 1.0e-12);
    CHECK(std::abs(sqe.values()[0] + sqe.values()[1]) < 1.0e-12);

    const auto sqeq0_transfer = ((-hydrogen_electronegativity + bond_interaction - source_h) -
                                 (-oxygen_electronegativity - bond_interaction - source_o)) /
                                denominator;
    const auto sqeq0 = calculate_pair(make_pair(1, -1), "sqeq0", 0.0, source);
    CHECK(std::abs(sqeq0.values()[0] - (1.0 + sqeq0_transfer)) < 1.0e-12);
    CHECK(std::abs(sqeq0.values()[1] - (-1.0 - sqeq0_transfer)) < 1.0e-12);
    CHECK(std::abs(sqeq0.values()[0] + sqeq0.values()[1]) < 1.0e-12);
}
