#include "bindings.h"
#include "native_parameter_catalog.h"

#include <chargefw/parameters/io/parameter_set_io.h>

#include <nanobind/stl/string.h>
#include <nanobind/stl/vector.h>

#include <filesystem>
#include <string>

namespace nb = nanobind;

namespace chargefw::python {
namespace {

auto make_parameter_catalog(const std::string& directory) -> NativeParameterCatalog {
    return NativeParameterCatalog{
        parameters::load_parameter_sets_json_directory(std::filesystem::path{directory})};
}

auto make_parameter_set(const std::string& path) -> NativeParameterCatalog {
    return NativeParameterCatalog{
        {parameters::load_parameter_set_json_file(std::filesystem::path{path})}};
}

auto make_parameter_sets(const std::string& directory) -> std::vector<NativeParameterCatalog> {
    auto result = std::vector<NativeParameterCatalog>{};
    for (auto& parameter_set :
         parameters::load_parameter_sets_json_directory(std::filesystem::path{directory})) {
        result.emplace_back(std::vector<parameters::ParameterSet>{std::move(parameter_set)});
    }
    return result;
}

auto combine_parameter_catalogs(const nb::list& catalogs) -> NativeParameterCatalog {
    auto parameter_sets = std::vector<parameters::ParameterSet>{};
    for (const auto& value : catalogs) {
        const auto& catalog = nb::cast<const NativeParameterCatalog&>(value);
        const auto& values = catalog.parameter_sets();
        parameter_sets.insert(parameter_sets.end(), values.begin(), values.end());
    }
    return NativeParameterCatalog{std::move(parameter_sets)};
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
        .def_prop_ro("size", &NativeParameterCatalog::size)
        .def("_descriptors", &parameter_set_descriptors);
    module.def("_load_parameter_catalog", &make_parameter_catalog, nb::arg("directory"));
    module.def("_load_parameter_set", &make_parameter_set, nb::arg("path"));
    module.def("_load_parameter_sets", &make_parameter_sets, nb::arg("directory"));
    module.def("_load_parameter_catalog_from_sets", &combine_parameter_catalogs,
               nb::arg("parameter_sets"));
}

} // namespace chargefw::python
