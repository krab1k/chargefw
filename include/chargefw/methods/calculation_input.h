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
                     const MethodOptions& method_options,
                     const features::ConformerFeatures* geometry = nullptr,
                     const parameters::ParameterView* parameters = nullptr)
        : prepared_molecule_{prepared_molecule}, method_options_{method_options},
          geometry_{geometry}, parameters_{parameters} {}

    CalculationInput(features::PreparedMolecule&&, const MethodOptions&,
                     const features::ConformerFeatures* = nullptr,
                     const parameters::ParameterView* = nullptr) = delete;

    CalculationInput(const features::PreparedMolecule&, MethodOptions&&,
                     const features::ConformerFeatures* = nullptr,
                     const parameters::ParameterView* = nullptr) = delete;

    [[nodiscard]] auto prepared_molecule() const noexcept -> const features::PreparedMolecule& {
        return prepared_molecule_.get();
    }

    [[nodiscard]] auto molecule() const noexcept -> const core::Molecule& {
        return prepared_molecule().molecule();
    }

    [[nodiscard]] auto topology() const noexcept -> const features::TopologyFeatures& {
        return prepared_molecule().topology();
    }

    [[nodiscard]] auto method_options() const noexcept -> const MethodOptions& {
        return method_options_.get();
    }

    [[nodiscard]] auto has_geometry() const noexcept -> bool {
        return geometry_ != nullptr;
    }

    [[nodiscard]] auto geometry_if_available() const noexcept
        -> const features::ConformerFeatures* {
        return geometry_;
    }

    [[nodiscard]] auto geometry() const -> const features::ConformerFeatures& {
        if (geometry_ == nullptr) {
            throw std::logic_error{"calculation input has no geometry"};
        }

        return *geometry_;
    }

    [[nodiscard]] auto has_parameters() const noexcept -> bool {
        return parameters_ != nullptr;
    }

    [[nodiscard]] auto parameters_if_available() const noexcept
        -> const parameters::ParameterView* {
        return parameters_;
    }

    [[nodiscard]] auto parameters() const -> const parameters::ParameterView& {
        if (parameters_ == nullptr) {
            throw std::logic_error{"calculation input has no parameters"};
        }

        return *parameters_;
    }

  private:
    std::reference_wrapper<const features::PreparedMolecule> prepared_molecule_;
    std::reference_wrapper<const MethodOptions> method_options_;
    const features::ConformerFeatures* geometry_ = nullptr;
    const parameters::ParameterView* parameters_ = nullptr;
};

} // namespace chargefw::methods