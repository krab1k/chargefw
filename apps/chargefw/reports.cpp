#include "cli_support.h"

#include <chargefw/core/periodic_table.h>
#include <chargefw/methods/method_registry.h>
#include <chargefw/parameters/io/parameter_set_io.h>

#include <algorithm>
#include <iostream>
#include <map>
#include <print>
#include <stdexcept>
#include <string_view>
#include <variant>
#include <vector>

namespace chargefw::cli {
namespace {

struct SummaryRow {
    std::string_view id;
    std::string_view name;
};

void print_summary(const std::vector<SummaryRow>& rows) {
    std::size_t id_width = 0;
    for (const auto& row : rows) {
        id_width = std::max(id_width, row.id.size());
    }
    for (const auto& row : rows) {
        std::println("{:<{}}  {}", row.id, id_width, row.name);
    }
}

} // namespace

void print_inspection(const ImportedCollection& imported) {
    std::println("records: {}", imported.molecules.size());
    for (std::size_t index = 0; index < imported.molecules.size(); ++index) {
        const auto& molecule = imported.molecules[index];
        std::map<int, std::size_t> elements;
        for (const auto& atom : molecule.atoms()) {
            ++elements[atom.atomic_number()];
        }
        std::print("record {} ({}) atoms={} bonds={} conformers={} coordinates={} formal_charge={} "
                   "elements=",
                   index, imported.export_context.records[index].identity.record_id,
                   molecule.atom_count(), molecule.bond_count(), molecule.conformer_count(),
                   molecule.has_coordinates(), core::total_formal_charge(molecule));
        bool first = true;
        for (const auto& [atomic_number, count] : elements) {
            std::print("{}{}:{}", first ? "" : ",", core::element_symbol(atomic_number), count);
            first = false;
        }
        std::println();
    }
}

void print_applicability(const calculation::AssessmentResult& assessment) {
    std::println("runnable plans: {}", assessment.plans().size());
    for (const auto& plan : assessment.plans()) {
        const auto& candidate = plan.candidate();
        std::println("plan method={} parameter_set={} execution={}", candidate.method->id(),
                     candidate.parameter_set == nullptr ? "-" : candidate.parameter_set->id(),
                     calculation::to_string(plan.policy().mode()));
    }
    std::println("rejected alternatives: {}", assessment.rejections().size());
    for (const auto& rejected : assessment.rejections()) {
        std::print("rejected method={}", rejected.method_id);
        if (rejected.parameter_set_id.has_value()) {
            std::print(" parameter_set={}", *rejected.parameter_set_id);
        }
        if (rejected.policy.has_value()) {
            std::print(" execution={}", calculation::to_string(rejected.policy->mode()));
        }
        for (const auto& issue : rejected.issues) {
            std::visit([](const auto& value) { std::print("; {}", value.message); }, issue);
        }
        std::println();
    }
    const auto* plan = assessment.default_plan();
    if (plan == nullptr) {
        std::println("selected execution: none");
        return;
    }
    const auto& selected = plan->candidate();
    std::println("selected method={} parameter_set={} execution={}", selected.method->id(),
                 selected.parameter_set == nullptr ? "-" : selected.parameter_set->id(),
                 calculation::to_string(plan->policy().mode()));
}

void print_methods(const std::string& method_id) {
    const auto& registry = methods::method_registry();
    if (method_id.empty()) {
        auto rows = std::vector<SummaryRow>{};
        rows.reserve(registry.methods().size());
        for (const auto& method : registry.methods()) {
            rows.push_back({method->id(), method->metadata().full_name});
        }
        print_summary(rows);
        return;
    }

    const auto* method = registry.find(method_id);
    if (method == nullptr) {
        throw std::invalid_argument{"method '" + method_id + "' is not registered"};
    }

    const auto& metadata = method->metadata();
    const auto requirements = method->requirements();
    std::println("id: {}", metadata.id);
    std::println("name: {}", metadata.name);
    std::println("full name: {}", metadata.full_name);
    std::println("publication: {}", metadata.publication.value_or("-"));
    std::println("notes: {}", metadata.notes.empty() ? "-" : metadata.notes);
    std::println("priority: {}", metadata.priority);
    std::println("requires coordinates: {}", requirements.coordinates ? "yes" : "no");
    std::println("time complexity: {}", methods::complexity_notation(requirements.resources.time));
    std::println("memory complexity: {}",
                 methods::complexity_notation(requirements.resources.memory));
    std::println("supports cutoff: {}", requirements.resources.supports_cutoff ? "yes" : "no");
    std::println("supports cover: {}", requirements.resources.supports_cover ? "yes" : "no");
    std::println("options:{}", method->option_schema().empty() ? " none" : "");
    for (const auto& option : method->option_schema()) {
        std::print("  {} (default=", option.id);
        std::visit([](const auto& value) { std::print("{}", value); }, option.default_value);
        std::print(")");
        if (!option.choices.empty()) {
            std::print(" choices=");
            for (std::size_t index = 0; index < option.choices.size(); ++index) {
                if (index != 0) {
                    std::print(",");
                }
                std::visit([](const auto& value) { std::print("{}", value); },
                           option.choices[index]);
            }
        }
        if (option.minimum.has_value()) {
            std::print("{}", option.minimum_inclusive ? " minimum>=" : " minimum>");
            std::visit([](const auto& value) { std::print("{}", value); }, *option.minimum);
        }
        if (option.maximum.has_value()) {
            std::print("{}", option.maximum_inclusive ? " maximum<=" : " maximum<");
            std::visit([](const auto& value) { std::print("{}", value); }, *option.maximum);
        }
        std::println();
    }
}

void print_parameter_sets(const std::string& parameter_set_id, const std::string& method_id) {
    const auto parameter_sets = parameters::load_default_parameter_sets();
    if (parameter_set_id.empty()) {
        if (!method_id.empty() && methods::method_registry().find(method_id) == nullptr) {
            throw std::invalid_argument{"method '" + method_id + "' is not registered"};
        }
        const auto visible = [&method_id](const parameters::ParameterSet& parameter_set) {
            return method_id.empty() || parameter_set.method_id() == method_id;
        };
        auto rows = std::vector<SummaryRow>{};
        rows.reserve(parameter_sets.size());
        for (const auto& parameter_set : parameter_sets) {
            if (visible(parameter_set)) {
                rows.push_back({parameter_set.id(), parameter_set.name()});
            }
        }
        print_summary(rows);
        return;
    }
    if (!method_id.empty()) {
        throw std::invalid_argument{"parameter-set ID and --method cannot be combined"};
    }

    const auto parameter_set =
        std::ranges::find(parameter_sets, parameter_set_id, &parameters::ParameterSet::id);
    if (parameter_set == parameter_sets.end()) {
        throw std::invalid_argument{"parameter set '" + parameter_set_id + "' was not found"};
    }

    const auto& metadata = parameter_set->metadata();
    std::println("id: {}", metadata.id);
    std::println("method: {}", metadata.method_id);
    std::println("name: {}", metadata.name);
    std::println("publication: {}", metadata.publication.empty() ? "-" : metadata.publication);
    std::println("notes: {}", metadata.notes.empty() ? "-" : metadata.notes);
    std::println("priority: {}", metadata.priority);
}

} // namespace chargefw::cli
