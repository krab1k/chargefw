#include <CLI/CLI.hpp>
#include <algorithm>
#include <chargefw/adapters/native/mol.h>
#include <chargefw/adapters/native/mol2.h>
#include <chargefw/adapters/native/sdf.h>
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
#include <string>
#include <utility>
#include <vector>

namespace charges = chargefw::charges;
namespace calculation = chargefw::calculation;
namespace core = chargefw::core;
namespace features = chargefw::features;
namespace methods = chargefw::methods;
namespace parameters = chargefw::parameters;
namespace adapters = chargefw::adapters;

namespace {

template <typename Reader>
auto read_collection(Reader& reader, const std::string& input_path) -> core::MoleculeCollection {
    std::vector<core::Molecule> molecules;

    while (const auto record = reader.next()) {
        if (!record->has_value()) {
            const auto& error = record->error();
            std::cerr << "Skipping record " << error.identity.record_index << ": " << error.message
                      << '\n';
            continue;
        }

        molecules.push_back(record->value().molecule);
    }

    if (molecules.empty()) {
        throw std::runtime_error{"No valid molecules found in input file: " + input_path};
    }

    return core::MoleculeCollection{std::move(molecules), input_path};
}

auto read_collection(const std::string& input_path) -> core::MoleculeCollection {
    std::ifstream input{input_path};
    if (!input) {
        throw std::runtime_error{"Unable to open input file: " + input_path};
    }

    auto extension = std::filesystem::path{input_path}.extension().string();
    std::ranges::transform(extension, extension.begin(), [](const unsigned char character) {
        return static_cast<char>(std::tolower(character));
    });
    if (extension == ".sdf") {
        adapters::native::sdf::SdfReader reader{input, input_path};
        return read_collection(reader, input_path);
    }
    if (extension == ".mol") {
        adapters::native::mol::MolReader reader{input, input_path};
        return read_collection(reader, input_path);
    }
    if (extension == ".mol2") {
        adapters::native::mol2::Mol2Reader reader{input, input_path};
        return read_collection(reader, input_path);
    }

    throw std::runtime_error{"Unsupported input file type: " + extension +
                             ". Supported types: .sdf, .mol, .mol2"};
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

auto print_charge_set(const charges::ChargeSet& charge_set) -> void {
    std::cout << "method: " << charge_set.method_id();

    if (charge_set.parameter_set_id().has_value()) {
        std::cout << "  parameters: " << *charge_set.parameter_set_id();
    }

    std::cout << '\n';

    for (const auto& assignment : charge_set.assignments()) {
        std::cout << "  molecule " << assignment.target.molecule_index;

        if (assignment.target.conformer_index.has_value()) {
            std::cout << " conformer " << *assignment.target.conformer_index;
        }

        std::cout << " charges:";

        for (const auto charge : assignment.charges.values()) {
            std::cout << ' ' << charge;
        }

        std::cout << "  total=" << assignment.charges.total() << '\n';
    }
}

} // namespace

auto main(int argc, char* argv[]) -> int {
    try {
        CLI::App app{"Calculate empirical partial atomic charges from a molecular file."};
        std::string input_path;
        app.add_option("input", input_path, "Input .sdf, .mol, or .mol2 file")->required();
        CLI11_PARSE(app, argc, argv);

        const auto collection = read_collection(input_path);
        const features::PreparedMoleculeCollection prepared_collection{collection};

        const auto parameter_sets = parameters::load_default_parameter_sets();

        const auto& registry = methods::method_registry();
        const auto candidates = method_pointers(registry);

        const auto result = calculation::calculate(
            calculation::CalculationRequest{.molecules = prepared_collection,
                                            .candidate_methods = candidates,
                                            .parameter_sets = parameter_sets});

        std::cout << "Loaded methods: " << candidates.size() << '\n';
        std::cout << "Loaded parameter sets: " << parameter_sets.size() << '\n';
        std::cout << "Loaded molecules: " << collection.size() << '\n';
        std::cout << "Applicable candidates: " << result.applicability.applicable.size() << '\n';

        if (!result.calculated()) {
            std::cerr << "No applicable method found\n";
            return 1;
        }

        if (const auto& charges = result.charges) {
            print_charge_set(*charges);
        }

        return 0;
    } catch (const std::exception& error) {
        std::cerr << "Fatal error: " << error.what() << '\n';
        return 1;
    }
}
