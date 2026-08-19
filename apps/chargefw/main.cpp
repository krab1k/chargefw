#include "cli_support.h"

#include <exception>
#include <iostream>
#include <print>
#include <stdexcept>

namespace {

auto run(int argc, char* argv[]) -> int {
    CLI::App app{"ChargeFW molecular charge calculation and inspection."};
    chargefw::cli::InputArguments calculate_input;
    chargefw::cli::InputArguments inspect_input;
    chargefw::cli::InputArguments applicability_input;
    chargefw::cli::SelectionArguments calculate_selection;
    chargefw::cli::SelectionArguments applicability_selection;
    std::string output_directory;
    std::string parameter_method;
    auto* calculate = app.add_subcommand("calculate", "Calculate and write partial charges");
    auto* inspect = app.add_subcommand("inspect", "Inspect imported molecular records");
    auto* applicability = app.add_subcommand("applicability", "Report applicable charge methods");
    auto* methods = app.add_subcommand("methods", "List registered charge methods");
    auto* parameters = app.add_subcommand("parameters", "List bundled parameter sets");
    chargefw::cli::add_input_options(*calculate, calculate_input);
    calculate->add_option("output", output_directory, "Output directory")->required();
    chargefw::cli::add_selection_options(*calculate, calculate_selection);
    chargefw::cli::add_input_options(*inspect, inspect_input);
    chargefw::cli::add_input_options(*applicability, applicability_input);
    chargefw::cli::add_selection_options(*applicability, applicability_selection);
    parameters->add_option("method", parameter_method, "Limit results to a method ID");
    CLI11_PARSE(app, argc, argv);

    if (*methods) {
        chargefw::cli::print_methods();
        return 0;
    }
    if (*parameters) {
        chargefw::cli::print_parameter_sets(parameter_method);
        return 0;
    }
    if (*inspect) {
        chargefw::cli::print_inspection(chargefw::cli::import_input(inspect_input));
        return 0;
    }
    if (*applicability) {
        const auto imported = chargefw::cli::import_input(applicability_input);
        chargefw::cli::print_applicability(chargefw::calculation::assess(
            chargefw::cli::make_request(imported, applicability_selection)));
        return 0;
    }
    if (!*calculate) {
        throw std::invalid_argument{
            "a subcommand is required; use calculate, inspect, applicability, "
            "methods, or parameters"};
    }

    const auto imported = chargefw::cli::import_input(calculate_input);
    const auto request = chargefw::cli::make_request(imported, calculate_selection);
    return chargefw::cli::write_calculation_outputs(output_directory, calculate_input.path,
                                                    imported, request,
                                                    chargefw::calculation::calculate(request));
}

} // namespace

auto main(int argc, char* argv[]) noexcept -> int {
    try {
        return run(argc, argv);
    } catch (const std::exception& error) {
        std::print(std::cerr, "Fatal error: {}\n", error.what());
        return 1;
    }
}
