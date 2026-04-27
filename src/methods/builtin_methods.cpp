#include "methods/builtin_methods.h"

#include "methods/builtin/dummy.h"
#include "methods/builtin/formal.h"

#include <memory>
#include <vector>

namespace chargefw::methods {

auto make_builtin_methods() -> std::vector<std::unique_ptr<Method>>
{
    std::vector<std::unique_ptr<Method>> methods;

    methods.push_back(std::make_unique<builtin::DummyMethod>());
    methods.push_back(std::make_unique<builtin::FormalMethod>());

    return methods;
}

} // namespace chargefw::methods