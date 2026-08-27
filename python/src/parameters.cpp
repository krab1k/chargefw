#include "bindings.h"
#include "native_parameter_catalog.h"

#include <chargefw/parameters/io/parameter_set_io.h>

#include <nanobind/stl/string.h>

#include <filesystem>
#include <string>

namespace nb = nanobind;

namespace chargefw::python {
namespace {

auto make_parameter_catalog(const std::string& directory) -> NativeParameterCatalog {
    return NativeParameterCatalog{
        parameters::load_parameter_sets_json_directory(std::filesystem::path{directory})};
}

} // namespace

void bind_parameters(nb::module_& module) {
    nb::class_<NativeParameterCatalog>(module, "_NativeParameterCatalog")
        .def_prop_ro("size", &NativeParameterCatalog::size);
    module.def("_load_parameter_catalog", &make_parameter_catalog, nb::arg("directory"));
}

} // namespace chargefw::python
