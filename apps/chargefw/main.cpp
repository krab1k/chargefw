#include "cli_support.h"

#include <algorithm>
#include <chrono>
#include <cstddef>
#include <exception>
#include <iostream>
#include <mutex>
#include <print>
#include <span>
#include <stdexcept>
#include <string>
#include <utility>

namespace {

class TerminalProgressObserver final : public chargefw::calculation::CalculationObserver {
  public:
    void on_progress(const chargefw::calculation::CalculationProgress& progress) const override {
        const std::scoped_lock lock{mutex_};
        switch (progress.phase) {
        case chargefw::calculation::CalculationPhase::computation_started:
            std::print(std::cerr, "Calculating {}...\n", progress.method_id);
            break;
        case chargefw::calculation::CalculationPhase::target_finished:
            ++completed_targets_;
            render("Targets", completed_targets_, progress.target_count);
            break;
        case chargefw::calculation::CalculationPhase::fragment_progress:
            render("Fragments", progress.completed_fragment_count, progress.fragment_count,
                   progress.target_index + 1, progress.target_count);
            break;
        case chargefw::calculation::CalculationPhase::computation_finished:
            std::print(std::cerr, "\n");
            break;
        default:
            break;
        }
        std::cerr.flush();
    }

  private:
    static constexpr std::size_t bar_width = 30;

    void render(const std::string_view label, const std::size_t complete, const std::size_t total,
                const std::optional<std::size_t> target_index = std::nullopt,
                const std::optional<std::size_t> target_count = std::nullopt) const {
        const auto filled = total == 0 ? 0U : std::min(bar_width, complete * bar_width / total);
        auto bar = std::string(bar_width, ' ');
        std::fill_n(bar.begin(), filled, '=');
        if (target_index.has_value() && target_count.has_value()) {
            std::print(std::cerr, "\r{} {}/{} [{}] {}/{}", label, *target_index, *target_count, bar,
                       complete, total);
            return;
        }
        std::print(std::cerr, "\r{} [{}] {}/{}", label, bar, complete, total);
    }

    mutable std::mutex mutex_;
    mutable std::size_t completed_targets_ = 0;
};

auto run(std::span<char*> arguments) -> int {
    CLI::App app{"ChargeFW molecular charge calculation and inspection."};
    chargefw::cli::InputArguments calculate_input;
    chargefw::cli::InputArguments inspect_input;
    chargefw::cli::InputArguments applicability_input;
    chargefw::cli::SelectionArguments calculate_selection;
    chargefw::cli::SelectionArguments applicability_selection;
    std::string output_directory;
    std::string parameter_method;
    bool progress = false;
    auto* calculate = app.add_subcommand("calculate", "Calculate and write partial charges");
    auto* inspect = app.add_subcommand("inspect", "Inspect imported molecular records");
    auto* applicability = app.add_subcommand("applicability", "Report applicable charge methods");
    auto* methods = app.add_subcommand("methods", "List registered charge methods");
    auto* parameters = app.add_subcommand("parameters", "List bundled parameter sets");
    chargefw::cli::add_input_options(*calculate, calculate_input);
    calculate->add_option("output", output_directory, "Output directory")->required();
    chargefw::cli::add_selection_options(*calculate, calculate_selection);
    calculate->add_flag("--progress", progress, "Show calculation progress on standard error");
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
        auto imported = chargefw::cli::import_input(applicability_input);
        chargefw::cli::print_applicability(chargefw::calculation::assess(
            chargefw::cli::make_request(std::move(imported.molecules), applicability_selection)));
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
    auto imported = chargefw::cli::import_input(calculate_input);
    run.metrics.parsing_seconds =
        std::chrono::duration<double>{std::chrono::steady_clock::now() - parsing_started}.count();
    auto export_context = std::move(imported.export_context);
    auto request = chargefw::cli::make_request(std::move(imported.molecules), calculate_selection);
    const auto requested_provenance =
        chargefw::cli::make_requested_provenance(export_context, request);
    const auto max_threads = request.resource_policy.max_threads;
    auto assessment = chargefw::calculation::assess(std::move(request));
    for (const auto& warning : assessment.execution_issues()) {
        std::println(std::cerr, "Warning: {}", warning.message);
    }
    const auto progress_observer = TerminalProgressObserver{};
    const auto& observer =
        progress ? static_cast<const chargefw::calculation::CalculationObserver&>(progress_observer)
                 : chargefw::calculation::default_calculation_observer();
    const auto result =
        chargefw::calculation::calculate(std::move(assessment), max_threads, observer);
    run.metrics.applicability_seconds = result.metrics.applicability_seconds;
    run.metrics.computation_seconds = result.metrics.computation_seconds;
    return chargefw::cli::write_calculation_outputs(
        output_directory, calculate_input.path, export_context, requested_provenance, result, run);
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
