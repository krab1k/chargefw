#pragma once

#include <string_view>
#include <optional>

namespace chargefw::methods {
struct MethodMetadata {
    std::string_view id;
    std::string_view name;
    std::string_view full_name;
    std::optional<std::string_view> publication;
    int priority = 0;
};

} // namespace chargefw::methods