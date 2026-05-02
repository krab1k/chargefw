#include <chargefw/parameters/models/parameter_key.h>

#include <stdexcept>
#include <string>

namespace chargefw::parameters {

auto atom_classification_kind_from_string(const std::string_view value)
    -> AtomParameterClassificationKind {
    if (value == "plain") {
        return AtomParameterClassificationKind::PLAIN;
    }

    if (value == "hbo") {
        return AtomParameterClassificationKind::HIGHEST_BOND_ORDER;
    }

    if (value == "bonded") {
        return AtomParameterClassificationKind::BONDED_ELEMENTS;
    }

    throw std::invalid_argument{"unknown atom parameter classification '" + std::string{value} +
                                "'"};
}

auto bond_classification_kind_from_string(const std::string_view value)
    -> BondParameterClassificationKind {
    if (value == "plain") {
        return BondParameterClassificationKind::PLAIN;
    }

    if (value == "bo") {
        return BondParameterClassificationKind::BOND_ORDER;
    }

    throw std::invalid_argument{"unknown bond parameter classification '" + std::string{value} +
                                "'"};
}

auto to_string(const AtomParameterClassificationKind kind) noexcept -> std::string_view {
    switch (kind) {
    case AtomParameterClassificationKind::PLAIN:
        return "plain";

    case AtomParameterClassificationKind::HIGHEST_BOND_ORDER:
        return "hbo";

    case AtomParameterClassificationKind::BONDED_ELEMENTS:
        return "bonded";
    }

    return "unknown";
}

auto to_string(const BondParameterClassificationKind kind) noexcept -> std::string_view {
    switch (kind) {
    case BondParameterClassificationKind::PLAIN:
        return "plain";

    case BondParameterClassificationKind::BOND_ORDER:
        return "bo";
    }

    return "unknown";
}

} // namespace chargefw::parameters