#include "bindings.h"

#include <chargefw/adapters/conformer_selection.h>
#include <chargefw/adapters/gemmi/input_options.h>
#include <chargefw/adapters/gemmi/mmcif_input.h>
#include <chargefw/adapters/gemmi/pdb_input.h>
#include <chargefw/adapters/molecule_record.h>
#include <chargefw/core/bond.h>

#include <nanobind/stl/array.h>
#include <nanobind/stl/string.h>
#include <nanobind/stl/vector.h>

#include <array>
#include <cstddef>
#include <sstream>
#include <string>
#include <utility>
#include <vector>

namespace nb = nanobind;

namespace chargefw::python {
namespace {

struct MoleculePayload {
    std::vector<int> atomic_numbers;
    std::vector<int> formal_charges;
    std::vector<std::array<std::size_t, 3>> bonds;
    std::vector<std::vector<std::array<double, 3>>> coordinates;
    std::string name;
    std::vector<std::string> atom_names;
    std::vector<std::string> conformer_names;
    adapters::MoleculeRecordIdentity identity;
    std::vector<adapters::MoleculeRecordDiagnostic> diagnostics;
};

auto make_payload(adapters::ImportedMoleculeRecord record) -> MoleculePayload {
    MoleculePayload result;
    const auto& molecule = record.molecule;
    result.atomic_numbers.reserve(molecule.atom_count());
    result.formal_charges.reserve(molecule.atom_count());
    result.atom_names.reserve(molecule.atom_count());
    for (const auto& atom : molecule.atoms()) {
        result.atomic_numbers.push_back(atom.atomic_number());
        result.formal_charges.push_back(atom.formal_charge());
        result.atom_names.emplace_back(atom.name());
    }

    result.bonds.reserve(molecule.bond_count());
    for (const auto& bond : molecule.bonds()) {
        result.bonds.push_back({bond.first_atom_index(), bond.second_atom_index(),
                                core::bond_order_value(bond.order())});
    }

    result.coordinates.reserve(molecule.conformer_count());
    result.conformer_names.reserve(molecule.conformer_count());
    for (const auto& conformer : molecule.conformers()) {
        auto positions = std::vector<std::array<double, 3>>{};
        positions.reserve(conformer.size());
        for (const auto& position : conformer.positions()) {
            positions.push_back({position.x, position.y, position.z});
        }
        result.coordinates.push_back(std::move(positions));
        result.conformer_names.emplace_back(conformer.name());
    }

    result.name = molecule.name();
    result.identity = std::move(record.identity);
    result.diagnostics = std::move(record.diagnostics);
    return result;
}

auto as_python(const MoleculePayload& payload) -> nb::dict {
    nb::dict result;
    result["atomic_numbers"] = nb::cast(payload.atomic_numbers);
    result["formal_charges"] = nb::cast(payload.formal_charges);
    result["bonds"] = nb::cast(payload.bonds);
    result["coordinates"] = nb::cast(payload.coordinates);
    result["name"] = payload.name;
    result["atom_names"] = nb::cast(payload.atom_names);
    result["conformer_names"] = nb::cast(payload.conformer_names);
    result["source"] = payload.identity.source;
    result["record_index"] = payload.identity.record_index;
    result["record_id"] = payload.identity.record_id;
    nb::list diagnostics;
    for (const auto& diagnostic : payload.diagnostics) {
        nb::dict value;
        value["code"] = diagnostic.code;
        value["message"] = diagnostic.message;
        value["line"] = diagnostic.line ? nb::cast(*diagnostic.line) : nb::none();
        diagnostics.append(std::move(value));
    }
    result["diagnostics"] = std::move(diagnostics);
    return result;
}

auto as_python(const std::vector<MoleculePayload>& payloads) -> nb::list {
    nb::list result;
    for (const auto& payload : payloads) {
        result.append(as_python(payload));
    }
    return result;
}

auto read_pdb(std::string contents, std::string source, const std::string& selection,
              const std::string& bonds, const std::string& conformers) -> nb::list {
    const auto native_selection = adapters::gemmi::record_selection_from_string(selection);
    const auto native_bonds = adapters::gemmi::bond_strategy_from_string(bonds);
    const auto native_conformers = adapters::conformer_selection_from_string(conformers);
    auto payloads = std::vector<MoleculePayload>{};
    {
        nb::gil_scoped_release release;
        std::istringstream input{std::move(contents)};
        auto reader = adapters::gemmi::pdb_input::PdbReader{input,
                                                            std::move(source),
                                                            {.selection = native_selection,
                                                             .bond_strategy = native_bonds,
                                                             .conformers = native_conformers}};
        if (auto record = reader.next()) {
            payloads.push_back(make_payload(std::move(*record)));
        }
    }
    return as_python(payloads);
}

auto read_mmcif(std::string contents, std::string source, const std::string& selection,
                const std::string& bonds, const std::string& conformers) -> nb::list {
    const auto native_selection = adapters::gemmi::record_selection_from_string(selection);
    const auto native_bonds = adapters::gemmi::bond_strategy_from_string(bonds);
    const auto native_conformers = adapters::conformer_selection_from_string(conformers);
    auto payloads = std::vector<MoleculePayload>{};
    {
        nb::gil_scoped_release release;
        std::istringstream input{std::move(contents)};
        auto reader = adapters::gemmi::mmcif_input::MmcifReader{input,
                                                                std::move(source),
                                                                {.selection = native_selection,
                                                                 .bond_strategy = native_bonds,
                                                                 .conformers = native_conformers}};
        while (auto record = reader.next()) {
            payloads.push_back(make_payload(std::move(*record)));
        }
    }
    return as_python(payloads);
}

} // namespace

void bind_adapters(nb::module_& module) {
    module.def("_read_pdb", &read_pdb, nb::arg("contents"), nb::arg("source"), nb::arg("selection"),
               nb::arg("bonds"), nb::arg("conformers"));
    module.def("_read_mmcif", &read_mmcif, nb::arg("contents"), nb::arg("source"),
               nb::arg("selection"), nb::arg("bonds"), nb::arg("conformers"));
}

} // namespace chargefw::python
