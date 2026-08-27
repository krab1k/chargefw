#include "bindings.h"

#include <chargefw/config.h>

#include <nanobind/nanobind.h>

namespace nb = nanobind;

NB_MODULE(_chargefw, module) {
    module.doc() = "Private native extension for ChargeFW.";
    module.def("version", [] { return CHARGEFW_VERSION_STRING; });

    auto core = module.def_submodule("core");
    auto methods = module.def_submodule("methods");
    auto parameters = module.def_submodule("parameters");
    auto calculation = module.def_submodule("calculation");

    chargefw::python::bind_core(core);
    chargefw::python::bind_methods(methods);
    chargefw::python::bind_parameters(parameters);
    chargefw::python::bind_calculation(calculation);
}
