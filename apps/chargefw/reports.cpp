#include "cli_support.h"

#include <chargefw/core/periodic_table.h>
#include <chargefw/methods/method_registry.h>
#include <chargefw/parameters/io/parameter_set_io.h>

#include <iostream>
#include <map>
#include <print>
#include <variant>

namespace chargefw::cli {

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
    const auto& report = assessment.applicability();
    std::println("applicable candidates: {}", report.applicable.size());
    for (const auto& candidate : report.applicable) {
        std::print("applicable method={} parameter_set={}", candidate.method_id,
                   candidate.parameter_set_id.value_or("-"));
        for (const auto& execution : candidate.execution_assessments) {
            std::print(" {}={}", calculation::to_string(execution.mode),
                       methods::to_string(execution.availability));
        }
        std::println();
    }
    std::println("rejected candidates: {}", report.rejected.size());
    for (const auto& rejected : report.rejected) {
        std::print("rejected method={}", rejected.method_id);
        if (rejected.parameter_set_id.has_value()) {
            std::print(" parameter_set={}", *rejected.parameter_set_id);
        }
        for (const auto& issue : rejected.issues) {
            std::print("; {}", issue.message);
        }
        std::println();
    }
    if (!report.selected_candidate_index.has_value()) {
        std::println("selected execution: none");
        return;
    }
    const auto& selected = report.applicable.at(*report.selected_candidate_index);
    std::println("selected method={} parameter_set={} execution={}", selected.method_id,
                 selected.parameter_set_id.value_or("-"),
                 calculation::to_string(assessment.execution_policy()->mode()));
}

void print_methods() {
    for (const auto& method : methods::method_registry().methods()) {
        std::println("{}\t{}", method->id(), method->metadata().name);
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
}

void print_parameter_sets(const std::string& method_id) {
    for (const auto& parameter_set : parameters::load_default_parameter_sets()) {
        if (method_id.empty() || parameter_set.method_id() == method_id) {
            std::println("{}\t{}\t{}", parameter_set.id(), parameter_set.method_id(),
                         parameter_set.name());
        }
    }
}

} // namespace chargefw::cli
