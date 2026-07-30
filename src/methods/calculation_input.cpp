#include <chargefw/core/molecule.h>
#include <chargefw/features/conformer_features.h>
#include <chargefw/features/prepared_molecule.h>
#include <chargefw/methods/calculation_input.h>
#include <chargefw/methods/method_options.h>

namespace chargefw::methods {

CalculationInput::CalculationInput(const features::PreparedMolecule& prepared_molecule,
                                   const MethodOptions& method_options, const double target_charge,
                                   const features::ConformerFeatures* geometry,
                                   const parameters::ParameterView* parameters)
    : prepared_molecule_{prepared_molecule}, method_options_{method_options},
      target_charge_{target_charge}, geometry_{geometry}, parameters_{parameters} {}
auto CalculationInput::prepared_molecule() const noexcept -> const features::PreparedMolecule& {
    return prepared_molecule_.get();
}

auto CalculationInput::molecule() const noexcept -> const core::Molecule& {
    return prepared_molecule().molecule();
}

auto CalculationInput::topology() const noexcept -> const features::TopologyFeatures& {
    return prepared_molecule().topology();
}

auto CalculationInput::method_options() const noexcept -> const MethodOptions& {
    return method_options_.get();
}

auto CalculationInput::target_charge() const noexcept -> double {
    return target_charge_;
}

auto CalculationInput::has_geometry() const noexcept -> bool {
    return geometry_ != nullptr;
}

auto CalculationInput::geometry_if_available() const noexcept
    -> const features::ConformerFeatures* {
    return geometry_;
}

auto CalculationInput::geometry() const -> const features::ConformerFeatures& {
    if (geometry_ == nullptr) {
        throw std::logic_error{"calculation input has no geometry"};
    }

    return *geometry_;
}

auto CalculationInput::has_parameters() const noexcept -> bool {
    return parameters_ != nullptr;
}

auto CalculationInput::parameters_if_available() const noexcept
    -> const parameters::ParameterView* {
    return parameters_;
}

auto CalculationInput::parameters() const -> const parameters::ParameterView& {
    if (parameters_ == nullptr) {
        throw std::logic_error{"calculation input has no parameters"};
    }

    return *parameters_;
}

} // namespace chargefw::methods
