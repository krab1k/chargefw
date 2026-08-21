#pragma once

#include <chargefw/features/prepared_molecule_collection.h>
#include <chargefw/methods/method_applicability.h>

#include <optional>
#include <string>

namespace chargefw::methods::detail {

auto validate_selected_candidate(const ApplicableMethod& selected,
                                 const features::PreparedMoleculeCollection& molecules) -> void;
auto validate_coordinate_targets(const ApplicableMethod& selected,
                                 const features::PreparedMoleculeCollection& molecules) -> void;
[[nodiscard]] auto parameter_set_id_for(const ApplicableMethod& selected)
    -> std::optional<std::string>;

} // namespace chargefw::methods::detail
