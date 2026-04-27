#include <chargefw/charges/charge_validation.h>

#include <stdexcept>
#include <string>

namespace chargefw::charges {
namespace {

auto make_context(const std::size_t charge_set_index, const std::size_t assignment_index)
    -> std::string {
    return "charge set " + std::to_string(charge_set_index) + ", assignment " +
           std::to_string(assignment_index) + ": ";
}

auto validate_target_molecule_index(const core::MoleculeCollection& molecules,
                                    const ChargeTarget& target, const std::size_t charge_set_index,
                                    const std::size_t assignment_index) -> void {
    if (target.molecule_index >= molecules.molecule_count()) {
        throw std::invalid_argument{make_context(charge_set_index, assignment_index) +
                                    "target molecule index is outside the molecule collection"};
    }
}

auto validate_target_conformer_index(const core::MoleculeCollection& molecules,
                                     const ChargeTarget& target, const std::size_t charge_set_index,
                                     const std::size_t assignment_index) -> void {
    if (!target.conformer_index.has_value()) {
        return;
    }

    const auto& molecule = molecules[target.molecule_index];

    if (*target.conformer_index >= molecule.conformer_count()) {
        throw std::invalid_argument{make_context(charge_set_index, assignment_index) +
                                    "target conformer index is outside the molecule"};
    }
}

auto validate_charge_count(const core::MoleculeCollection& molecules,
                           const ChargeAssignment& assignment, const std::size_t charge_set_index,
                           const std::size_t assignment_index) -> void {
    const auto& molecule = molecules[assignment.target.molecule_index];

    if (assignment.charges.size() != molecule.atom_count()) {
        throw std::invalid_argument{make_context(charge_set_index, assignment_index) +
                                    "atomic charge count does not match molecule atom count"};
    }
}

auto validate_assignment(const core::MoleculeCollection& molecules,
                         const ChargeAssignment& assignment, const std::size_t charge_set_index,
                         const std::size_t assignment_index) -> void {
    validate_target_molecule_index(molecules, assignment.target, charge_set_index,
                                   assignment_index);
    validate_target_conformer_index(molecules, assignment.target, charge_set_index,
                                    assignment_index);
    validate_charge_count(molecules, assignment, charge_set_index, assignment_index);
}

} // namespace

auto validate_charge_collection(const core::MoleculeCollection& molecules,
                                const ChargeCollection& charges) -> void {
    const auto charge_sets = charges.charge_sets();

    for (std::size_t charge_set_index = 0; charge_set_index < charge_sets.size();
         ++charge_set_index) {
        const auto assignments = charge_sets[charge_set_index].assignments();

        for (std::size_t assignment_index = 0; assignment_index < assignments.size();
             ++assignment_index) {
            validate_assignment(molecules, assignments[assignment_index], charge_set_index,
                                assignment_index);
        }
    }
}

} // namespace chargefw::charges