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

[[nodiscard]] inline auto single_bond_key(const int first_atomic_number,
                                          const int second_atomic_number)
    -> parameters::BondParameterKey {
    return {.first_atom = plain_atom_key(first_atomic_number),
            .second_atom = plain_atom_key(second_atomic_number),
            .bond = {.classification = parameters::BondParameterClassificationKind::BOND_ORDER,
                     .type = "1"}};
}

} // namespace chargefw::test
