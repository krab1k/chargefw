#pragma once

#include <filesystem>
#include <vector>

namespace chargefw::parameters {

[[nodiscard]] auto bundled_parameter_directory() -> std::filesystem::path;

[[nodiscard]] auto default_parameter_directories() -> std::vector<std::filesystem::path>;

} // namespace chargefw::parameters