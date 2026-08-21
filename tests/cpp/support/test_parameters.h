#pragma once

#include <chargefw/parameters/models/parameter_key.h>

#include <string>
#include <utility>

namespace chargefw::test {

[[nodiscard]] inline auto atom_key(const int atomic_number,
                                   const parameters::AtomParameterClassificationKind classification,
                                   std::string type) -> parameters::AtomParameterKey {
    return {
        .atomic_number = atomic_number, .classification = classification, .type = std::move(type)};
}

[[nodiscard]] inline auto plain_atom_key(const int atomic_number) -> parameters::AtomParameterKey {
    return atom_key(atomic_number, parameters::AtomParameterClassificationKind::PLAIN, "*");
}

[[nodiscard]] inline auto hbo_atom_key(const int atomic_number, std::string type)
    -> parameters::AtomParameterKey {
    return atom_key(atomic_number, parameters::AtomParameterClassificationKind::HIGHEST_BOND_ORDER,
                    std::move(type));
}

[[nodiscard]] inline auto single_bond_key(const int first_atomic_number,
                                          const int second_atomic_number)
    -> parameters::BondParameterKey {
    return {.first_atom = plain_atom_key(first_atomic_number),
            .second_atom = plain_atom_key(second_atomic_number),
            .bond = {.classification = parameters::BondParameterClassificationKind::BOND_ORDER,
                     .type = "1"}};
}

/// Undirected plain bond key: matches a bond between the two elements regardless of direction.
[[nodiscard]] inline auto plain_bond_key(const int first_atomic_number,
                                         const int second_atomic_number)
    -> parameters::BondParameterKey {
    return {.first_atom = plain_atom_key(first_atomic_number),
            .second_atom = plain_atom_key(second_atomic_number),
            .bond = {.classification = parameters::BondParameterClassificationKind::PLAIN,
                     .type = "*"}};
}

/// Wildcard plain bond key: matches any bond regardless of element or direction.
[[nodiscard]] inline auto plain_bond_key() -> parameters::BondParameterKey {
    return {.first_atom = plain_atom_key(0),
            .second_atom = plain_atom_key(0),
            .bond = {.classification = parameters::BondParameterClassificationKind::PLAIN,
                     .type = "*"}};
}

} // namespace chargefw::test
