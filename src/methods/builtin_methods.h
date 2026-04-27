#pragma once

#include <chargefw/methods/method.h>

#include <memory>
#include <vector>

namespace chargefw::methods {

[[nodiscard]] auto make_builtin_methods() -> std::vector<std::unique_ptr<Method>>;

} // namespace chargefw::methods