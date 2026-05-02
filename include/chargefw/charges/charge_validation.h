#pragma once

#include <chargefw/charges/charge_collection.h>
#include <chargefw/core/molecule_collection.h>

namespace chargefw::charges {

auto validate_charge_collection(const core::MoleculeCollection& collection,
                                const ChargeCollection& charges) -> void;

} // namespace chargefw::charges