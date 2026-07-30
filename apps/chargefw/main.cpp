#include <CLI/CLI.hpp>
#include <algorithm>
#include <chargefw/adapters/native/json_input.h>
#include <chargefw/adapters/native/json_output.h>
#include <chargefw/adapters/native/mol2_input.h>
#include <chargefw/adapters/native/mol_input.h>
#include <chargefw/adapters/native/sdf_input.h>
#include <chargefw/calculation/calculation.h>
#include <chargefw/charges/charge_collection.h>
#include <chargefw/core/molecule.h>
#include <chargefw/core/molecule_collection.h>
#include <chargefw/features/prepared_molecule_collection.h>
#include <chargefw/methods/method.h>
#include <chargefw/methods/method_registry.h>
#include <chargefw/parameters/io/parameter_set_io.h>

#include <cctype>
#include <exception>
#include <filesystem>
#include <fstream>
#include <iostream>
#include <optional>
#include <print>
#include <stdexcept>
#include <string>
#include <utility>
#include <vector>

namespace calculation = chargefw::calculation;
namespace core = chargefw::core;
namespace features = chargefw::features;
namespace methods = chargefw::methods;
namespace parameters = chargefw::parameters;
namespace adapters = chargefw::adapters;

namespace {

struct ImportedCollection {
    core::MoleculeCollection molecules;
    std::vector<adapters::MoleculeRecordIdentity> identities;
    std::vector<adapters::MoleculeRecordMapping> mappings;
};

template <typename Reader>
auto read_collection(Reader& reader, const std::string& input_path) -> ImportedCollection {
    std::vector<core::Molecule> molecules;
    std::vector<adapters::MoleculeRecordIdentity> identities;
    std::vector<adapters::MoleculeRecordMapping> mappings;

    while (const auto record = reader.next()) {
        auto imported = std::move(*record);
        molecules.push_back(std::move(imported.molecule));
        identities.push_back(std::move(imported.identity));
        mappings.push_back(std::move(imported.mapping));
    }

    if (molecules.empty()) {
        throw std::runtime_error{"No molecules found in input file: " + input_path};
    }

    return ImportedCollection{.molecules =
                                  core::MoleculeCollection{std::move(molecules), input_path},
                              .identities = std::move(identities),
                              .mappings = std::move(mappings)};
}

auto read_collection(const std::string& input_path) -> ImportedCollection {
    std::ifstream input{input_path};
    if (!input) {
        throw std::runtime_error{"Unable to open input file: " + input_path};
    }

    auto extension = std::filesystem::path{input_path}.extension().string();
    std::ranges::transform(extension, extension.begin(), [](const unsigned char character) -> char {
        return static_cast<char>(std::tolower(character));
    });
    if (extension == ".sdf") {
        adapters::native::sdf_input::SdfReader reader{input, input_path};
        return read_collection(reader, input_path);
    }
    if (extension == ".mol") {
        adapters::native::mol_input::MolReader reader{input, input_path};
        return read_collection(reader, input_path);
    }
    if (extension == ".mol2") {
        adapters::native::mol2_input::Mol2Reader reader{input, input_path};
        return read_collection(reader, input_path);
    }
    if (extension == ".json") {
        adapters::native::json_input::JsonReader reader{input, input_path};
        return read_collection(reader, input_path);
    }

    throw std::runtime_error{"Unsupported input file type: " + extension +
                             ". Supported types: .sdf, .mol, .mol2, .json"};
}

auto method_pointers(const methods::MethodRegistry& registry)
    -> std::vector<const methods::Method*> {
    std::vector<const methods::Method*> result;
    result.reserve(registry.methods().size());

    for (const auto& method : registry.methods()) {
        result.push_back(method.get());
    }

    return result;
}

[[nodiscard]] auto result_document(const ImportedCollection& imported,
                                   const calculation::CalculationResult& result)
    -> adapters::ChargeResultDocument {
    auto document = adapters::ChargeResultDocument{
        .generator_name = "ChargeFW", .generator_version = "0.0.1", .records = {}};
    document.records.reserve(imported.molecules.size());
    for (std::size_t index = 0; index < imported.molecules.size(); ++index) {
        document.records.push_back(
            adapters::ChargeResultRecord{.identity = imported.identities[index],
                                         .mapping = imported.mappings[index],
                                         .charges = result.charges});
    }
    return document;
}

} // namespace

auto main(int argc, char* argv[]) -> int {
    try {
        CLI::App app{"Calculate empirical partial atomic charges from a molecular file."};
        std::string input_path;
        std::string output_path;
        app.add_option("input", input_path, "Input .sdf, .mol, .mol2, or ChargeFW .json file")
            ->required();
        app.add_option("-o,--output", output_path, "Write JSON calculation result to this file");
        CLI11_PARSE(app, argc, argv);

        const auto imported = read_collection(input_path);
        const features::PreparedMoleculeCollection prepared_collection{imported.molecules};

        const auto parameter_sets = parameters::load_default_parameter_sets();

        const auto& registry = methods::method_registry();
        const auto candidates = method_pointers(registry);

        const auto result = calculation::calculate(
            calculation::CalculationRequest{.molecules = prepared_collection,
                                            .candidate_methods = candidates,
                                            .parameter_sets = parameter_sets});

        const auto output = result_document(imported, result);
        if (output_path.empty()) {
            adapters::native::json_output::JsonWriter{std::cout}.write(output);
        } else {
            std::ofstream output_file{output_path};
            if (!output_file) {
                throw std::runtime_error{"Unable to open output file: " + output_path};
            }
            adapters::native::json_output::JsonWriter{output_file}.write(output);
        }

        if (!result.calculated()) {
            return 1;
        }

        return 0;
    } catch (const std::exception& error) {
        std::print(std::cerr, "Fatal error: {}\n", error.what());
        return 1;
    }
}
