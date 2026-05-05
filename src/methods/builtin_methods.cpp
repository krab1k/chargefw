#include "methods/builtin_methods.h"

#include "methods/builtin/charge2.h"
#include "methods/builtin/delre.h"
#include "methods/builtin/dummy.h"
#include "methods/builtin/formal.h"
#include "methods/builtin/gdac.h"
#include "methods/builtin/mgc.h"
#include "methods/builtin/mpeoe.h"
#include "methods/builtin/peoe.h"
#include "methods/builtin/veem.h"

#include <memory>
#include <vector>

namespace chargefw::methods {

auto make_builtin_methods() -> std::vector<std::unique_ptr<Method>> {
    std::vector<std::unique_ptr<Method>> methods;

    methods.push_back(std::make_unique<builtin::DummyMethod>());
    methods.push_back(std::make_unique<builtin::FormalMethod>());
    methods.push_back(std::make_unique<builtin::VEEMMethod>());
    methods.push_back(std::make_unique<builtin::PEOEMethod>());
    methods.push_back(std::make_unique<builtin::MPEOEMethod>());
    methods.push_back(std::make_unique<builtin::GDACMethod>());
    methods.push_back(std::make_unique<builtin::Charge2Method>());
    methods.push_back(std::make_unique<builtin::DelReMethod>());
    methods.push_back(std::make_unique<builtin::MGCMethod>());

    return methods;
}

} // namespace chargefw::methods