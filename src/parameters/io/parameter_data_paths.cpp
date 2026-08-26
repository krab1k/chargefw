#include <chargefw/config.h>
#include <chargefw/parameters/io/parameter_data_paths.h>

#include <dlfcn.h>
#include <filesystem>
#include <stdexcept>
#include <vector>

namespace chargefw::parameters {
namespace {

[[nodiscard]] auto loaded_library_directory() -> std::filesystem::path {
    Dl_info information{};
    if (dladdr(reinterpret_cast<const void*>(&bundled_parameter_directory), &information) == 0 ||
        information.dli_fname == nullptr) {
        throw std::runtime_error{"could not determine the ChargeFW library location"};
    }

    return std::filesystem::absolute(information.dli_fname).parent_path();
}

} // namespace

auto bundled_parameter_directory() -> std::filesystem::path {
    return (loaded_library_directory() / CHARGEFW_PARAMETER_DIRECTORY_FROM_LIBRARY)
        .lexically_normal();
}

auto default_parameter_directories() -> std::vector<std::filesystem::path> {
    return {bundled_parameter_directory()};
}

} // namespace chargefw::parameters
