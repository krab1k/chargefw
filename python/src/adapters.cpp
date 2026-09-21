#include "bindings.h"
#include "native_execution_result.h"
#include "native_input_metadata.h"

#include <chargefw/adapters/charge_result.h>
#include <chargefw/adapters/conformer_selection.h>
#include <chargefw/adapters/gemmi/input_options.h>
#include <chargefw/adapters/gemmi/mmcif_input.h>
#include <chargefw/adapters/gemmi/mmcif_output.h>
#include <chargefw/adapters/gemmi/pdb_input.h>
#include <chargefw/adapters/generated_output.h>
#include <chargefw/adapters/molecule_record.h>
#include <chargefw/adapters/native/json_input.h>
#include <chargefw/adapters/native/json_output.h>
#include <chargefw/adapters/native/mol2_input.h>
#include <chargefw/adapters/native/mol_input.h>
#include <chargefw/adapters/native/sdf_input.h>
#include <chargefw/calculation/calculation.h>
#include <chargefw/config.h>
#include <chargefw/core/bond.h>
#include <chargefw/core/molecule.h>

#include <nanobind/stl/array.h>
#include <nanobind/stl/optional.h>
#include <nanobind/stl/shared_ptr.h>
#include <nanobind/stl/string.h>
#include <nanobind/stl/tuple.h>
#include <nanobind/stl/vector.h>

#include <array>
#include <cstddef>
#include <memory>
#include <optional>
#include <sstream>
#include <stdexcept>
#include <string>
#include <string_view>
#include <tuple>
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
    std::optional<adapters::MoleculeImportMetadata> import_metadata;
    std::shared_ptr<NativeInputMetadata> native_input_metadata;
};

[[nodiscard]] auto format_name(const adapters::MolecularSourceFormat format) -> std::string_view {
    switch (format) {
    case adapters::MolecularSourceFormat::molecule_json:
        return "molecule-json";
    case adapters::MolecularSourceFormat::mol:
        return "mol";
    case adapters::MolecularSourceFormat::sdf:
        return "sdf";
    case adapters::MolecularSourceFormat::mol2:
        return "mol2";
    case adapters::MolecularSourceFormat::pdb:
        return "pdb";
    case adapters::MolecularSourceFormat::mmcif:
        return "mmcif";
    }
    throw std::invalid_argument{"unknown molecular source format"};
}

[[nodiscard]] auto connectivity_name(const adapters::SourceConnectivity connectivity)
    -> std::string_view {
    switch (connectivity) {
    case adapters::SourceConnectivity::absent:
        return "absent";
    case adapters::SourceConnectivity::explicitly_empty:
        return "explicitly-empty";
    case adapters::SourceConnectivity::present:
        return "present";
    }
    throw std::invalid_argument{"unknown source connectivity state"};
}

[[nodiscard]] auto atom_references(const std::vector<adapters::SourceAtomReference>& references)
    -> nb::list {
    auto result = nb::list{};
    for (const auto& reference : references) {
        if (!reference.structural_labels.has_value()) {
            result.append(nb::make_tuple(reference.position, reference.id, nb::none()));
            continue;
        }
        const auto hierarchy = [](const adapters::SourceHierarchyLabels& labels) {
            auto value = nb::dict{};
            value["atom"] = labels.atom;
            value["residue"] = labels.residue;
            value["chain"] = labels.chain;
            value["sequence"] = labels.sequence;
            return value;
        };
        const auto& labels = *reference.structural_labels;
        auto encoded = nb::dict{};
        encoded["author"] = hierarchy(labels.author);
        encoded["label"] = hierarchy(labels.label);
        encoded["entity"] = labels.entity;
        encoded["insertion_code"] = labels.insertion_code;
        encoded["alternate_location"] = labels.alternate_location;
        encoded["segment"] = labels.segment;
        result.append(nb::make_tuple(reference.position, reference.id, std::move(encoded)));
    }
    return result;
}

[[nodiscard]] auto import_metadata(const adapters::MoleculeImportMetadata& metadata) -> nb::dict {
    auto conformers = nb::list{};
    for (const auto& conformer : metadata.conformers) {
        conformers.append(
            nb::make_tuple(conformer.position, conformer.id, atom_references(conformer.sites)));
    }
    auto result = nb::dict{};
    result["format"] = std::string{format_name(metadata.format)};
    result["atoms"] = atom_references(metadata.atoms);
    result["conformers"] = std::move(conformers);
    result["record_selection"] = metadata.record_selection;
    result["alternate_location_selection"] = metadata.alternate_location_selection;
    result["conformer_selection"] = metadata.conformer_selection;
    result["bond_strategy"] = metadata.bond_strategy;
    result["source_connectivity"] = std::string{connectivity_name(metadata.source_connectivity)};
    return result;
}

[[nodiscard]] auto portable_id(const adapters::PortableId& id) -> nb::object {
    if (id.empty()) {
        return nb::none();
    }
    return std::visit([](const auto& value) { return nb::cast(value); }, *id.value());
}

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
    result.native_input_metadata = std::make_shared<NativeInputMetadata>(
        record.identity, record.diagnostics, record.import_metadata);
    result.identity = std::move(record.identity);
    result.diagnostics = std::move(record.diagnostics);
    result.import_metadata = std::move(record.import_metadata);
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
    result["record_id"] = portable_id(payload.identity.record_id);
    auto diagnostics = nb::list{};
    for (const auto& diagnostic : payload.diagnostics) {
        diagnostics.append(nb::make_tuple(diagnostic.code, diagnostic.message, diagnostic.line));
    }
    result["diagnostics"] = std::move(diagnostics);
    result["native_input_metadata"] = payload.native_input_metadata;
    if (payload.import_metadata.has_value()) {
        result["import_metadata"] = import_metadata(*payload.import_metadata);
    } else {
        result["import_metadata"] = nb::none();
    }
    return result;
}

auto as_python(const std::vector<MoleculePayload>& payloads) -> nb::list {
    nb::list result;
    for (const auto& payload : payloads) {
        result.append(as_python(payload));
    }
    return result;
}

template <typename Reader> auto read_all(Reader& reader) -> std::vector<MoleculePayload> {
    auto payloads = std::vector<MoleculePayload>{};
    while (auto record = reader.next()) {
        payloads.push_back(make_payload(std::move(*record)));
    }
    return payloads;
}

auto parse(std::string contents, std::string source, const std::string& format,
           const std::string& selection, const std::string& bonds, const std::string& conformers)
    -> nb::list {
    auto payloads = std::vector<MoleculePayload>{};
    {
        nb::gil_scoped_release release;
        std::istringstream input{std::move(contents)};
        if (format == "mol") {
            auto reader = adapters::native::mol_input::MolReader{input, std::move(source)};
            payloads = read_all(reader);
        } else if (format == "sdf") {
            auto reader = adapters::native::sdf_input::SdfReader{input, std::move(source)};
            payloads = read_all(reader);
        } else if (format == "mol2") {
            auto reader = adapters::native::mol2_input::Mol2Reader{input, std::move(source)};
            payloads = read_all(reader);
        } else if (format == "molecule-json") {
            auto reader = adapters::native::json_input::JsonReader{
                input, std::move(source), adapters::conformer_selection_from_string(conformers)};
            payloads = read_all(reader);
        } else if (format == "pdb" || format == "mmcif") {
            const auto options = adapters::gemmi::InputOptions{
                .selection = adapters::gemmi::record_selection_from_string(selection),
                .bond_strategy = adapters::gemmi::bond_strategy_from_string(bonds),
                .conformers = adapters::conformer_selection_from_string(conformers)};
            if (format == "pdb") {
                auto reader =
                    adapters::gemmi::pdb_input::PdbReader{input, std::move(source), options};
                payloads = read_all(reader);
            } else {
                auto reader =
                    adapters::gemmi::mmcif_input::MmcifReader{input, std::move(source), options};
                payloads = read_all(reader);
            }
        } else {
            throw std::invalid_argument{"unsupported molecular input format: " + format};
        }
    }
    return as_python(payloads);
}

auto dumps(const NativeExecutionResult& native_result, const std::string& format,
           const std::string& sdf_version) -> std::string {
    const auto& result = native_result.result();
    auto output = std::ostringstream{};
    {
        nb::gil_scoped_release release;
        if (format == "result-json") {
            adapters::native::json_output::JsonWriter{output}.write(result, "ChargeFW",
                                                                    CHARGEFW_VERSION_STRING);
        } else {
            if (!result.execution().calculated()) {
                throw std::invalid_argument{"molecular output requires a successful calculation"};
            }
            const auto output_format = [&format, &sdf_version] {
                if (format == "sdf") {
                    return sdf_version == "v2000" ? adapters::generated_output::Format::sdf_v2000
                                                  : adapters::generated_output::Format::sdf_v3000;
                }
                if (format == "mol2") {
                    return adapters::generated_output::Format::mol2;
                }
                if (format == "mmcif") {
                    return adapters::generated_output::Format::mmcif;
                }
                throw std::invalid_argument{"unsupported calculation output format: " + format};
            }();
            adapters::generated_output::write(output, result.inputs(), *result.execution().charges,
                                              output_format, "ChargeFW", CHARGEFW_VERSION_STRING);
        }
    }
    return output.str();
}

auto attach_mmcif(std::string contents, const NativeExecutionResult& native_result,
                  const nb::sequence& molecules, const std::string& selection,
                  const std::string& conformers, const bool overwrite) -> std::string {
    auto source_molecules = std::vector<core::Molecule>{};
    source_molecules.reserve(static_cast<std::size_t>(nb::len(molecules)));
    for (const auto molecule : molecules) {
        source_molecules.push_back(nb::cast<const core::Molecule&>(molecule));
    }
    const auto native_selection = adapters::gemmi::record_selection_from_string(selection);
    const auto native_conformers = adapters::conformer_selection_from_string(conformers);
    const auto& result = native_result.result().execution();
    auto output = std::ostringstream{};
    {
        nb::gil_scoped_release release;
        if (!result.calculated()) {
            throw std::invalid_argument{"charge attachment requires a successful calculation"};
        }
        auto input = std::istringstream{std::move(contents)};
        auto reader = adapters::gemmi::mmcif_input::MmcifReader{
            input,
            {},
            {.selection = native_selection,
             .bond_strategy = adapters::gemmi::BondStrategy::none,
             .conformers = native_conformers}};
        auto records = std::vector<adapters::ImportedMoleculeRecord>{};
        while (auto record = reader.next()) {
            records.push_back(std::move(*record));
        }
        if (records.size() != source_molecules.size()) {
            throw std::invalid_argument{
                "Gemmi target molecule count does not match the calculation input"};
        }
        for (std::size_t molecule_index = 0; molecule_index < records.size(); ++molecule_index) {
            const auto source_atoms = source_molecules[molecule_index].atoms();
            const auto target_atoms = records[molecule_index].molecule.atoms();
            if (source_atoms.size() != target_atoms.size()) {
                throw std::invalid_argument{
                    "Gemmi target atom count does not match the calculation input"};
            }
            for (std::size_t atom_index = 0; atom_index < source_atoms.size(); ++atom_index) {
                if (source_atoms[atom_index].atomic_number() !=
                        target_atoms[atom_index].atomic_number() ||
                    source_atoms[atom_index].formal_charge() !=
                        target_atoms[atom_index].formal_charge()) {
                    throw std::invalid_argument{
                        "Gemmi target atom order does not match the calculation input"};
                }
            }
        }
        const auto source_document = reader.source_document();
        const auto& block_indices = reader.source_block_indices();
        if (!overwrite) {
            for (const auto block_index : block_indices) {
                const auto& block = source_document->blocks.at(block_index);
                if (block.has_mmcif_category("_sb_ncbr_partial_atomic_charges_meta.") ||
                    block.has_mmcif_category("_sb_ncbr_partial_atomic_charges.")) {
                    throw std::invalid_argument{
                        "Gemmi target already contains partial charge categories"};
                }
            }
        }
        const auto source =
            adapters::gemmi::mmcif_output::MmcifSource{.document = source_document,
                                                       .block_indices = block_indices,
                                                       .selection = native_selection};
        adapters::gemmi::mmcif_output::MmcifWriter{output}.write_mmcif(
            records, *result.charges, source, "ChargeFW", CHARGEFW_VERSION_STRING,
            adapters::gemmi::mmcif_output::WriteMode::replace);
    }
    return output.str();
}

} // namespace

void bind_adapters(nb::module_& module) {
    [[maybe_unused]] const auto native_input_metadata =
        nb::class_<NativeInputMetadata>(module, "_NativeInputMetadata");
    module.def("_parse", &parse, nb::arg("contents"), nb::arg("source"), nb::arg("format"),
               nb::arg("selection"), nb::arg("bonds"), nb::arg("conformers"));
    module.def("_dumps", &dumps, nb::arg("result"), nb::arg("format"), nb::arg("sdf_version"));
    module.def("_attach_mmcif", &attach_mmcif, nb::arg("contents"), nb::arg("result"),
               nb::arg("molecules"), nb::arg("selection"), nb::arg("conformers"),
               nb::arg("overwrite"));
}

} // namespace chargefw::python
