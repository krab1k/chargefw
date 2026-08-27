#include <chargefw/config.h>

#include <nanobind/nanobind.h>

namespace nb = nanobind;

NB_MODULE(_chargefw, module) {
    module.doc() = "Private native extension for ChargeFW.";
    module.def("version", [] { return CHARGEFW_VERSION_STRING; });
}
