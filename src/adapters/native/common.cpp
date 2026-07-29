#include "common.h"

#include <charconv>
#include <sstream>
#include <stdexcept>
#include <string>
#include <utility>

namespace chargefw::adapters::native::common {

auto read_line(std::istream& input, std::size_t& line, const std::string_view record_name)
    -> std::string {
    std::string result;

    if (!std::getline(input, result)) {
        throw std::runtime_error{"unexpected end of " + std::string{record_name} + " record"};
    }

    ++line;
    return result;
}

auto trim(const std::string_view value) -> std::string_view {
    const auto first = value.find_first_not_of(" \t\r");

    if (first == std::string_view::npos) {
        return {};
    }

    return value.substr(first, value.find_last_not_of(" \t\r") - first + 1);
}

auto parse_int(const std::string_view value, const std::string_view field) -> int {
    int result = 0;
    const auto [end, error] = std::from_chars(value.begin(), value.end(), result);

    if (error != std::errc{} || end != value.end()) {
        throw std::runtime_error{"invalid " + std::string{field}};
    }

    return result;
}

auto parse_trimmed_int(const std::string_view value, const std::string_view field) -> int {
    return parse_int(trim(value), field);
}

auto parse_double(const std::string_view value, const std::string_view field) -> double {
    std::istringstream input{std::string{value}};
    double result = 0.0;
    char remaining = '\0';

    if (!(input >> result) || (input >> remaining)) {
        throw std::runtime_error{"invalid " + std::string{field}};
    }

    return result;
}

auto parse_trimmed_double(const std::string_view value, const std::string_view field) -> double {
    return parse_double(trim(value), field);
}

auto fixed_field(const std::string_view line, const std::size_t offset, const std::size_t width,
                 const std::string_view field) -> std::string_view {
    if (line.size() < offset + width) {
        throw std::runtime_error{"missing " + std::string{field}};
    }

    return line.substr(offset, width);
}

auto numeric_bond_order(const int value) -> core::BondOrder {
    switch (value) {
    case 1:
        return core::BondOrder::SINGLE;
    case 2:
        return core::BondOrder::DOUBLE;
    case 3:
        return core::BondOrder::TRIPLE;
    case 4:
        return core::BondOrder::AROMATIC;
    default:
        throw std::runtime_error{"unsupported bond order " + std::to_string(value)};
    }
}

auto identity_mapping(const std::size_t count) -> std::vector<std::optional<std::size_t>> {
    std::vector<std::optional<std::size_t>> mapping;
    mapping.reserve(count);

    for (std::size_t index = 0; index < count; ++index) {
        mapping.emplace_back(index);
    }

    return mapping;
}

} // namespace chargefw::adapters::native::common
