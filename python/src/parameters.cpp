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
    nb::gil_scoped_release release;
    return NativeParameterCatalog{
        parameters::load_parameter_sets_json_directory(std::filesystem::path{directory})};
}

auto parameter_set_descriptors(const NativeParameterCatalog& catalog) -> nb::list {
    auto result = nb::list{};
    for (const auto& metadata : catalog.descriptors()) {
        auto descriptor = nb::dict{};
        descriptor["id"] = metadata.id;
        descriptor["method_id"] = metadata.method_id;
        descriptor["name"] = metadata.name;
        descriptor["publication"] = metadata.publication;
        descriptor["notes"] = metadata.notes;
        descriptor["priority"] = metadata.priority;
        result.append(std::move(descriptor));
    }
    return result;
}

} // namespace

void bind_parameters(nb::module_& module) {
    nb::class_<NativeParameterCatalog>(module, "_NativeParameterCatalog")
        .def("_descriptors", &parameter_set_descriptors);
    module.def("_load_parameter_catalog", &make_parameter_catalog, nb::arg("directory"));
}

} // namespace chargefw::python
