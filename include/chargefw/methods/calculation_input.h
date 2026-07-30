#pragma once

#include <chargefw/features/conformer_features.h>
#include <chargefw/features/prepared_molecule.h>
#include <chargefw/methods/method_options.h>
#include <chargefw/parameters/models/parameter_view.h>

#include <functional>
#include <stdexcept>

namespace chargefw::methods {

class CalculationInput {
  public:
    CalculationInput(const features::PreparedMolecule& prepared_molecule,
                     const MethodOptions& method_options, double target_charge,
                     const features::ConformerFeatures* geometry = nullptr,
                     const parameters::ParameterView* parameters = nullptr);

    CalculationInput(features::PreparedMolecule&&, const MethodOptions&, double,
                     const features::ConformerFeatures* = nullptr,
                     const parameters::ParameterView* = nullptr) = delete;

    CalculationInput(const features::PreparedMolecule&, MethodOptions&&, double,
                     const features::ConformerFeatures* = nullptr,
                     const parameters::ParameterView* = nullptr) = delete;

    [[nodiscard]] auto prepared_molecule() const noexcept -> const features::PreparedMolecule&;

    [[nodiscard]] auto molecule() const noexcept -> const core::Molecule&;

    [[nodiscard]] auto topology() const noexcept -> const features::TopologyFeatures&;

    [[nodiscard]] auto method_options() const noexcept -> const MethodOptions&;

    [[nodiscard]] auto target_charge() const noexcept -> double;

    [[nodiscard]] auto has_geometry() const noexcept -> bool;

    [[nodiscard]] auto geometry_if_available() const noexcept -> const features::ConformerFeatures*;

    [[nodiscard]] auto geometry() const -> const features::ConformerFeatures&;

    [[nodiscard]] auto has_parameters() const noexcept -> bool;

    [[nodiscard]] auto parameters_if_available() const noexcept -> const parameters::ParameterView*;

    [[nodiscard]] auto parameters() const -> const parameters::ParameterView&;

  private:
    std::reference_wrapper<const features::PreparedMolecule> prepared_molecule_;
    std::reference_wrapper<const MethodOptions> method_options_;
    double target_charge_;
    const features::ConformerFeatures* geometry_ = nullptr;
    const parameters::ParameterView* parameters_ = nullptr;
};

} // namespace chargefw::methods
