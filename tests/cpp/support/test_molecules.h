#pragma once

#include <chargefw/core/molecule.h>

namespace chargefw::test {

[[nodiscard]] auto make_water() -> core::Molecule;

[[nodiscard]] auto make_formally_charged_pair() -> core::Molecule;

} // namespace chargefw::test