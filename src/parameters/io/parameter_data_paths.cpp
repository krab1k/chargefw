#include <chargefw/config.h>
#include <chargefw/parameters/io/parameter_data_paths.h>

#include <filesystem>
#include <vector>

namespace chargefw::parameters {
auto bundled_parameter_directory() -> std::filesystem::path {
    return std::filesystem::path{CHARGEFW_INSTALL_PARAMETER_DIR};
}

auto default_parameter_directories() -> std::vector<std::filesystem::path> {
    const auto installed_directory = bundled_parameter_directory();
    if (std::filesystem::exists(installed_directory)) {
        return {installed_directory};
    }

    return {std::filesystem::path{CHARGEFW_BUILD_PARAMETER_DIR}};
}

} // namespace chargefw::parameters
