#include <chargefw/parameters/io/parameter_data_paths.h>
#include <chargefw/config.h>

#include <cstdlib>
#include <filesystem>
#include <vector>

namespace chargefw::parameters {
namespace {

[[nodiscard]] auto environment_path(const char* name) -> std::filesystem::path {
    const auto* value = std::getenv(name);

    if (value == nullptr || *value == '\0') {
        return {};
    }

    return std::filesystem::path{value};
}

} // namespace

auto bundled_parameter_directory() -> std::filesystem::path {
    return std::filesystem::path{CHARGEFW_INSTALL_PARAMETER_DIR};
}

auto default_parameter_directories() -> std::vector<std::filesystem::path> {
    std::vector<std::filesystem::path> directories;

    const auto override_directory = environment_path("CHARGEFW_PARAMETER_DIR");

    if (!override_directory.empty()) {
        directories.push_back(override_directory);
    }

    directories.push_back(bundled_parameter_directory());

    return directories;
}

} // namespace chargefw::parameters