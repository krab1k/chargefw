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
                   index, imported.records[index].identity.record_id, molecule.atom_count(),
                   molecule.bond_count(), molecule.conformer_count(), molecule.has_coordinates(),
                   core::total_formal_charge(molecule));
        bool first = true;
        for (const auto& [atomic_number, count] : elements) {
            std::print("{}{}:{}", first ? "" : ",", core::element_symbol(atomic_number), count);
            first = false;
        }
        std::println();
    }
}

void print_applicability(const calculation::AssessmentResult& assessment) {
    std::println("applicable candidates: {}", assessment.applicability.applicable.size());
    for (const auto& candidate : assessment.applicability.applicable) {
        std::print("applicable method={} parameter_set={}", candidate.method->id(),
                   candidate.parameter_set == nullptr ? "-" : candidate.parameter_set->id());
        for (const auto& execution : candidate.execution_assessments) {
            const auto mode = execution.mode == calculation::ExecutionMode::full     ? "full"
                              : execution.mode == calculation::ExecutionMode::cutoff ? "cutoff"
                                                                                     : "cover";
            std::print(" {}={}", mode,
                       execution.availability == methods::ExecutionAvailability::unsupported
                           ? "unsupported"
                           : "available");
        }
        std::println();
    }
    std::println("rejected candidates: {}", assessment.applicability.rejected.size());
    const auto& registry = methods::method_registry();
    for (const auto& rejected : assessment.applicability.rejected) {
        std::print("rejected method={}", registry.methods()[rejected.method_index]->id());
        if (rejected.parameter_set_index.has_value()) {
            std::print(" parameter_set={}",
                       assessment.parameter_sets[*rejected.parameter_set_index].id());
        }
        for (const auto& issue : rejected.issues) {
            std::print("; {}", issue.message);
        }
        std::println();
    }
    if (assessment.selected == nullptr) {
        std::println("selected execution: none");
        return;
    }
    std::println(
        "selected method={} parameter_set={} execution={}", assessment.selected->method->id(),
        assessment.selected->parameter_set == nullptr ? "-"
                                                      : assessment.selected->parameter_set->id(),
        assessment.execution_policy->mode() == calculation::ExecutionMode::full     ? "full"
        : assessment.execution_policy->mode() == calculation::ExecutionMode::cutoff ? "cutoff"
                                                                                    : "cover");
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
