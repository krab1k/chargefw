#pragma once

#include <chargefw/methods/method.h>

#include <memory>
#include <span>
#include <string>
#include <string_view>
#include <vector>

namespace chargefw::methods {

class MethodRegistry {
public:
    explicit MethodRegistry(std::vector<std::unique_ptr<Method>> methods);

    [[nodiscard]] auto find(std::string_view id) const noexcept -> const Method*;

    [[nodiscard]] auto methods() const noexcept -> std::span<const std::unique_ptr<Method>>;

    [[nodiscard]] auto names() const -> std::vector<std::string>;

private:
    std::vector<std::unique_ptr<Method>> methods_;
};

[[nodiscard]] auto method_registry() -> const MethodRegistry&;

} // namespace chargefw::methods