#pragma once

#include <chargefw/parameters/parameter_set.h>

#include <filesystem>
#include <istream>
#include <vector>

namespace chargefw::parameters {

[[nodiscard]] auto load_parameter_set_json(std::istream& input) -> ParameterSet;

[[nodiscard]] auto load_parameter_set_json_file(const std::filesystem::path& path) -> ParameterSet;

[[nodiscard]] auto load_parameter_sets_json_directory(const std::filesystem::path& directory)
    -> std::vector<ParameterSet>;

[[nodiscard]] auto load_default_parameter_sets() -> std::vector<ParameterSet>;

} // namespace chargefw::parameters