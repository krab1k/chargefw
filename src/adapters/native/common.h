#pragma once

#include <chargefw/adapters/molecule_record.h>
#include <chargefw/core/bond.h>

#include <cstddef>
#include <istream>
#include <optional>
#include <string>
#include <string_view>
#include <vector>

namespace chargefw::adapters::native::common {

[[nodiscard]] auto read_line(std::istream& input, std::size_t& line, std::string_view record_name)
    -> std::string;

[[nodiscard]] auto trim(std::string_view value) -> std::string_view;

[[nodiscard]] auto parse_int(std::string_view value, std::string_view field) -> int;

[[nodiscard]] auto parse_trimmed_int(std::string_view value, std::string_view field) -> int;

[[nodiscard]] auto parse_double(std::string_view value, std::string_view field) -> double;

[[nodiscard]] auto parse_trimmed_double(std::string_view value, std::string_view field) -> double;

[[nodiscard]] auto fixed_field(std::string_view line, std::size_t offset, std::size_t width,
                               std::string_view field) -> std::string_view;

[[nodiscard]] auto numeric_bond_order(int value) -> core::BondOrder;

[[nodiscard]] auto identity_mapping(std::size_t count) -> std::vector<std::optional<std::size_t>>;

} // namespace chargefw::adapters::native::common
