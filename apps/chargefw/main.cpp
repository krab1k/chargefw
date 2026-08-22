#include "cli_support.h"

#include <chrono>
#include <cstddef>
#include <exception>
#include <iostream>
#include <print>
#include <span>
#include <stdexcept>
#include <utility>

namespace {

auto run(std::span<char*> arguments) -> int {
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
    const auto argc = static_cast<int>(arguments.size());
    auto* const argv = arguments.data();
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

    auto run =
        chargefw::cli::CalculationRun{.metrics = {}, .started = std::chrono::steady_clock::now()};
    run.metrics.started_at = chargefw::cli::utc_timestamp();
    const auto parsing_started = std::chrono::steady_clock::now();
    const auto imported = chargefw::cli::import_input(calculate_input);
    run.metrics.parsing_seconds =
        std::chrono::duration<double>{std::chrono::steady_clock::now() - parsing_started}.count();
    const auto request = chargefw::cli::make_request(imported, calculate_selection);
    auto assessment = chargefw::calculation::assess(request);
    for (const auto& warning : assessment.execution_issues) {
        std::println(std::cerr, "Warning: {}", warning.message);
    }
    const auto result = chargefw::calculation::calculate(std::move(assessment),
                                                         request.resource_policy.max_threads);
    run.metrics.applicability_seconds = result.metrics.applicability_seconds;
    run.metrics.computation_seconds = result.metrics.computation_seconds;
    return chargefw::cli::write_calculation_outputs(output_directory, calculate_input.path,
                                                    imported, request, result, run);
}

} // namespace

auto main(int argc, char* argv[]) noexcept -> int {
    try {
        return run({argv, static_cast<std::size_t>(argc)});
    } catch (const std::exception& error) {
        std::print(std::cerr, "Fatal error: {}\n", error.what());
        return 1;
    }
}
