#include <chargefw/parameters/parameter_set_io.h>

#include <chargefw/core/periodic_table.h>
#include <chargefw/parameters/parameter_data_paths.h>

#include <nlohmann/json.hpp>

#include <algorithm>
#include <cctype>
#include <cmath>
#include <filesystem>
#include <fstream>
#include <memory>
#include <set>
#include <stdexcept>
#include <string>
#include <string_view>
#include <utility>
#include <vector>

namespace chargefw::parameters {
namespace {

using Json = nlohmann::json;

[[nodiscard]] auto json_type_name(const Json& value) -> std::string {
    return std::string{value.type_name()};
}

[[noreturn]] auto throw_error(const std::string& context, const std::string& message) -> void {
    if (context.empty()) {
        throw std::invalid_argument{message};
    }

    throw std::invalid_argument{context + ": " + message};
}

[[nodiscard]] auto child_context(const std::string& context, const std::string_view child)
    -> std::string {
    if (context.empty()) {
        return std::string{child};
    }

    return context + "." + std::string{child};
}

[[nodiscard]] auto array_context(const std::string& context, const std::size_t index)
    -> std::string {
    return context + "[" + std::to_string(index) + "]";
}

auto ensure_object(const Json& value, const std::string& context) -> void {
    if (!value.is_object()) {
        throw_error(context, "expected object, got " + json_type_name(value));
    }
}

auto ensure_array(const Json& value, const std::string& context) -> void {
    if (!value.is_array()) {
        throw_error(context, "expected array, got " + json_type_name(value));
    }
}

[[nodiscard]] auto require_object(const Json& value, const std::string& context) -> const Json& {
    ensure_object(value, context);
    return value;
}

[[nodiscard]] auto require_array(const Json& value, const std::string& context) -> const Json& {
    ensure_array(value, context);
    return value;
}

[[nodiscard]] auto required_member(const Json& object, const std::string_view name,
                                   const std::string& context) -> const Json& {
    ensure_object(object, context);

    const auto found = object.find(std::string{name});

    if (found == object.end()) {
        throw_error(context, "missing required member '" + std::string{name} + "'");
    }

    return *found;
}

[[nodiscard]] auto optional_member(const Json& object, const std::string_view name) -> const Json* {
    if (!object.is_object()) {
        return nullptr;
    }

    const auto found = object.find(std::string{name});

    if (found == object.end()) {
        return nullptr;
    }

    return std::addressof(*found);
}

[[nodiscard]] auto require_string(const Json& value, const std::string& context) -> std::string {
    if (!value.is_string()) {
        throw_error(context, "expected string, got " + json_type_name(value));
    }

    return value.get<std::string>();
}

[[nodiscard]] auto require_double(const Json& value, const std::string& context) -> double {
    if (!value.is_number()) {
        throw_error(context, "expected number, got " + json_type_name(value));
    }

    const auto number = value.get<double>();

    if (!std::isfinite(number)) {
        throw_error(context, "expected finite number");
    }

    return number;
}

[[nodiscard]] auto optional_string_member(const Json& object, const std::string_view name,
                                          const std::string& context) -> std::string {
    const auto* value = optional_member(object, name);

    if (value == nullptr) {
        return {};
    }

    return require_string(*value, child_context(context, name));
}

[[nodiscard]] auto slugify(std::string value) -> std::string {
    std::string result;
    result.reserve(value.size());

    auto previous_was_separator = false;

    for (const auto character : value) {
        const auto unsigned_character = static_cast<unsigned char>(character);

        if (std::isalnum(unsigned_character) != 0) {
            result.push_back(static_cast<char>(std::tolower(unsigned_character)));
            previous_was_separator = false;
            continue;
        }

        if (!previous_was_separator && !result.empty()) {
            result.push_back('-');
            previous_was_separator = true;
        }
    }

    while (!result.empty() && result.back() == '-') {
        result.pop_back();
    }

    return result;
}

[[nodiscard]] auto generated_id_from_metadata(const std::string& method_id, const std::string& name)
    -> std::string {
    auto id = slugify(method_id + "-" + name);

    if (!id.empty()) {
        return id;
    }

    return slugify(method_id);
}

[[nodiscard]] auto parse_metadata(const Json& root, const std::string& fallback_id,
                                  const std::string& context) -> ParameterSetMetadata {
    const auto metadata_context = child_context(context, "metadata");
    const auto& metadata_json = required_member(root, "metadata", context);

    ensure_object(metadata_json, metadata_context);

    auto metadata = ParameterSetMetadata{};

    metadata.id = optional_string_member(metadata_json, "id", metadata_context);

    if (metadata.id.empty() && !fallback_id.empty()) {
        metadata.id = fallback_id;
    }

    if (const auto* method_id = optional_member(metadata_json, "method_id")) {
        metadata.method_id =
            require_string(*method_id, child_context(metadata_context, "method_id"));
    } else {
        metadata.method_id =
            require_string(required_member(metadata_json, "method", metadata_context),
                           child_context(metadata_context, "method"));
    }

    metadata.name = optional_string_member(metadata_json, "name", metadata_context);
    metadata.publication = optional_string_member(metadata_json, "publication", metadata_context);
    metadata.notes = optional_string_member(metadata_json, "notes", metadata_context);

    if (metadata.id.empty()) {
        metadata.id = generated_id_from_metadata(metadata.method_id, metadata.name);
    }

    return metadata;
}

[[nodiscard]] auto parse_names(const Json& object, const std::string& context)
    -> std::vector<std::string> {
    const auto names_context = child_context(context, "names");
    const auto& names_member = required_member(object, "names", context);
    const auto& names_json = require_array(names_member, names_context);

    std::vector<std::string> names;
    names.reserve(names_json.size());

    for (std::size_t index = 0; index < names_json.size(); ++index) {
        names.push_back(require_string(names_json[index], array_context(names_context, index)));
    }

    return names;
}

[[nodiscard]] auto make_named_parameters(const std::vector<std::string>& names,
                                         const Json& values_json, const std::string& context)
    -> std::vector<NamedParameter> {
    ensure_array(values_json, context);

    if (values_json.size() != names.size()) {
        throw_error(context, "value count " + std::to_string(values_json.size()) +
                                 " does not match name count " + std::to_string(names.size()));
    }

    std::vector<NamedParameter> parameters;
    parameters.reserve(names.size());

    for (std::size_t index = 0; index < names.size(); ++index) {
        parameters.push_back(NamedParameter{
            .name = names[index],
            .value = require_double(values_json[index], array_context(context, index))});
    }

    return parameters;
}

[[nodiscard]] auto parse_common_parameters(const Json& root, const std::string& context)
    -> CommonParameters {
    const auto* common_json = optional_member(root, "common");

    if (common_json == nullptr) {
        return CommonParameters{};
    }

    const auto common_context = child_context(context, "common");
    ensure_object(*common_json, common_context);

    const auto names = parse_names(*common_json, common_context);

    const auto values_context = child_context(common_context, "values");
    const auto& values_json = required_member(*common_json, "values", common_context);

    return CommonParameters{make_named_parameters(names, values_json, values_context)};
}
[[nodiscard]] auto atomic_number_from_symbol(const std::string& symbol, const std::string& context)
    -> int {
    if (symbol == "*") {
        return 0;
    }

    try {
        return core::periodic_table().element(symbol).atomic_number;
    } catch (const std::out_of_range&) {
        throw_error(context, "unknown element symbol '" + symbol + "'");
    }
}

[[nodiscard]] auto parse_atom_key(const Json& key_json, const std::string& context)
    -> AtomParameterKey {
    ensure_array(key_json, context);

    if (key_json.size() != 3) {
        throw_error(context, "atom key must have exactly 3 items: [symbol, classification, type]");
    }

    const auto symbol = require_string(key_json[0], array_context(context, 0));
    const auto classification = require_string(key_json[1], array_context(context, 1));
    const auto type = require_string(key_json[2], array_context(context, 2));

    return AtomParameterKey{.atomic_number =
                                atomic_number_from_symbol(symbol, array_context(context, 0)),
                            .classification = atom_classification_kind_from_string(classification),
                            .type = type};
}

[[nodiscard]] auto parse_bond_type_key(const Json& key_json, const std::string& context)
    -> BondTypeKey {
    ensure_array(key_json, context);

    if (key_json.size() != 2) {
        throw_error(context, "bond type key must have exactly 2 items: [classification, type]");
    }

    const auto classification = require_string(key_json[0], array_context(context, 0));
    const auto type = require_string(key_json[1], array_context(context, 1));

    return BondTypeKey{.classification = bond_classification_kind_from_string(classification),
                       .type = type};
}

[[nodiscard]] auto parse_flat_bond_key(const Json& key_json, const std::string& context)
    -> BondParameterKey {
    ensure_array(key_json, context);

    if (key_json.size() != 8) {
        throw_error(context, "flat bond key must have exactly 8 items: "
                             "[first_symbol, first_classification, first_type, "
                             "second_symbol, second_classification, second_type, "
                             "bond_classification, bond_type]");
    }

    const auto first_symbol = require_string(key_json[0], array_context(context, 0));
    const auto first_classification = require_string(key_json[1], array_context(context, 1));
    const auto first_type = require_string(key_json[2], array_context(context, 2));

    const auto second_symbol = require_string(key_json[3], array_context(context, 3));
    const auto second_classification = require_string(key_json[4], array_context(context, 4));
    const auto second_type = require_string(key_json[5], array_context(context, 5));

    const auto bond_classification = require_string(key_json[6], array_context(context, 6));
    const auto bond_type = require_string(key_json[7], array_context(context, 7));

    return BondParameterKey{
        .first_atom =
            AtomParameterKey{
                .atomic_number = atomic_number_from_symbol(first_symbol, array_context(context, 0)),
                .classification = atom_classification_kind_from_string(first_classification),
                .type = first_type},
        .second_atom = AtomParameterKey{.atomic_number = atomic_number_from_symbol(
                                            second_symbol, array_context(context, 3)),
                                        .classification = atom_classification_kind_from_string(
                                            second_classification),
                                        .type = second_type},
        .bond =
            BondTypeKey{.classification = bond_classification_kind_from_string(bond_classification),
                        .type = bond_type}};
}

[[nodiscard]] auto parse_bond_key(const Json& key_json, const std::string& context)
    -> BondParameterKey {
    ensure_array(key_json, context);

    if (key_json.size() == 8 && key_json[0].is_string()) {
        return parse_flat_bond_key(key_json, context);
    }

    if (key_json.size() != 3) {
        throw_error(context, "bond key must either have nested form "
                             "[first_atom, second_atom, bond] or flat ChargeFW2 form with 8 items");
    }

    return BondParameterKey{.first_atom = parse_atom_key(key_json[0], array_context(context, 0)),
                            .second_atom = parse_atom_key(key_json[1], array_context(context, 1)),
                            .bond = parse_bond_type_key(key_json[2], array_context(context, 2))};
}

[[nodiscard]] auto parse_atom_parameters(const Json& root, const std::string& context)
    -> AtomParameters {
    const auto* atom_json = optional_member(root, "atom");

    if (atom_json == nullptr) {
        return AtomParameters{};
    }

    const auto atom_context = child_context(context, "atom");
    ensure_object(*atom_json, atom_context);

    const auto names = parse_names(*atom_json, atom_context);

    const auto data_context = child_context(atom_context, "data");
    const auto& data_member = required_member(*atom_json, "data", atom_context);
    const auto& data_json = require_array(data_member, data_context);

    std::vector<AtomParameterEntry> entries;
    entries.reserve(data_json.size());

    for (std::size_t entry_index = 0; entry_index < data_json.size(); ++entry_index) {
        const auto entry_context = array_context(data_context, entry_index);
        const auto& entry_json = require_object(data_json[entry_index], entry_context);

        const auto key_context = child_context(entry_context, "key");
        const auto value_context = child_context(entry_context, "value");

        const auto key =
            parse_atom_key(required_member(entry_json, "key", entry_context), key_context);

        auto parameters = make_named_parameters(
            names, required_member(entry_json, "value", entry_context), value_context);

        entries.push_back(AtomParameterEntry{.key = key, .parameters = std::move(parameters)});
    }

    return AtomParameters{std::move(entries)};
}

[[nodiscard]] auto parse_bond_parameters(const Json& root, const std::string& context)
    -> BondParameters {
    const auto* bond_json = optional_member(root, "bond");

    if (bond_json == nullptr) {
        return BondParameters{};
    }

    const auto bond_context = child_context(context, "bond");
    ensure_object(*bond_json, bond_context);

    const auto names = parse_names(*bond_json, bond_context);

    const auto data_context = child_context(bond_context, "data");
    const auto& data_member = required_member(*bond_json, "data", bond_context);
    const auto& data_json = require_array(data_member, data_context);

    std::vector<BondParameterEntry> entries;
    entries.reserve(data_json.size());

    for (std::size_t entry_index = 0; entry_index < data_json.size(); ++entry_index) {
        const auto entry_context = array_context(data_context, entry_index);
        const auto& entry_json = require_object(data_json[entry_index], entry_context);

        const auto key_context = child_context(entry_context, "key");
        const auto value_context = child_context(entry_context, "value");

        const auto key =
            parse_bond_key(required_member(entry_json, "key", entry_context), key_context);

        auto parameters = make_named_parameters(
            names, required_member(entry_json, "value", entry_context), value_context);

        entries.push_back(BondParameterEntry{.key = key, .parameters = std::move(parameters)});
    }

    return BondParameters{std::move(entries)};
}

[[nodiscard]] auto parse_parameter_set(const Json& root, const std::string& fallback_id,
                                       const std::string& context) -> ParameterSet {
    ensure_object(root, context);

    return ParameterSet{parse_metadata(root, fallback_id, context),
                        parse_common_parameters(root, context),
                        parse_atom_parameters(root, context), parse_bond_parameters(root, context)};
}

[[nodiscard]] auto load_parameter_set_json_document(std::istream& input,
                                                    const std::string& fallback_id,
                                                    const std::string& context) -> ParameterSet {
    try {
        const auto root = Json::parse(input);
        return parse_parameter_set(root, fallback_id, context);
    } catch (const nlohmann::json::exception& error) {
        throw_error(context, error.what());
    }
}

[[nodiscard]] auto json_files_in_directory(const std::filesystem::path& directory)
    -> std::vector<std::filesystem::path> {
    if (!std::filesystem::exists(directory)) {
        throw std::invalid_argument{"parameter directory '" + directory.string() +
                                    "' does not exist"};
    }

    if (!std::filesystem::is_directory(directory)) {
        throw std::invalid_argument{"parameter path '" + directory.string() +
                                    "' is not a directory"};
    }

    std::vector<std::filesystem::path> files;

    for (const auto& entry : std::filesystem::directory_iterator{directory}) {
        if (!entry.is_regular_file()) {
            continue;
        }

        if (entry.path().extension() == ".json") {
            files.push_back(entry.path());
        }
    }

    std::ranges::sort(files);

    return files;
}

auto reject_duplicate_parameter_set_id(const std::set<std::string>& ids,
                                       const ParameterSet& parameter_set,
                                       const std::filesystem::path& source) -> void {
    if (ids.contains(std::string{parameter_set.id()})) {
        throw std::invalid_argument{"duplicate parameter set id '" +
                                    std::string{parameter_set.id()} + "' while loading '" +
                                    source.string() + "'"};
    }
}

auto append_parameter_sets_from_directory_if_present(std::vector<ParameterSet>& parameter_sets,
                                                     std::set<std::string>& ids,
                                                     const std::filesystem::path& directory)
    -> void {
    if (!std::filesystem::exists(directory)) {
        return;
    }

    if (!std::filesystem::is_directory(directory)) {
        throw std::invalid_argument{"default parameter path '" + directory.string() +
                                    "' is not a directory"};
    }

    for (const auto& file : json_files_in_directory(directory)) {
        auto parameter_set = load_parameter_set_json_file(file);
        reject_duplicate_parameter_set_id(ids, parameter_set, file);
        ids.insert(std::string{parameter_set.id()});
        parameter_sets.push_back(std::move(parameter_set));
    }
}

} // namespace

auto load_parameter_set_json(std::istream& input) -> ParameterSet {
    return load_parameter_set_json_document(input, {}, "parameter JSON");
}

auto load_parameter_set_json_file(const std::filesystem::path& path) -> ParameterSet {
    std::ifstream input{path};

    if (!input) {
        throw std::invalid_argument{"failed to open parameter file '" + path.string() + "'"};
    }

    return load_parameter_set_json_document(input, path.stem().string(),
                                            "parameter file '" + path.string() + "'");
}

auto load_parameter_sets_json_directory(const std::filesystem::path& directory)
    -> std::vector<ParameterSet> {
    std::vector<ParameterSet> parameter_sets;
    std::set<std::string> ids;

    for (const auto& file : json_files_in_directory(directory)) {
        auto parameter_set = load_parameter_set_json_file(file);
        reject_duplicate_parameter_set_id(ids, parameter_set, file);
        ids.insert(std::string{parameter_set.id()});
        parameter_sets.push_back(std::move(parameter_set));
    }

    return parameter_sets;
}

auto load_default_parameter_sets() -> std::vector<ParameterSet> {
    std::vector<ParameterSet> parameter_sets;
    std::set<std::string> ids;

    for (const auto& directory : default_parameter_directories()) {
        append_parameter_sets_from_directory_if_present(parameter_sets, ids, directory);
    }

    return parameter_sets;
}

} // namespace chargefw::parameters