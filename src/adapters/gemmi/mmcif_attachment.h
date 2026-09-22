#pragma once

#include <chargefw/adapters/charge_result.h>

#include <iosfwd>
#include <string_view>

#include <gemmi/cifdoc.hpp>

namespace chargefw::adapters::gemmi::mmcif_output {

auto write_attached(std::ostream& output, const ChargeCalculationResult& result,
                    const ::gemmi::cif::Document& source, bool overwrite,
                    std::string_view generator_version) -> void;

} // namespace chargefw::adapters::gemmi::mmcif_output
