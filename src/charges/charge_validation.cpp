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

auto validate_target_molecule_index(const core::MoleculeCollection& collection,
                                    const ChargeTarget& target, const std::size_t charge_set_index,
                                    const std::size_t assignment_index) -> void {
    if (target.molecule_index >= collection.size()) {
        throw std::invalid_argument{make_context(charge_set_index, assignment_index) +
                                    "target molecule index is outside the molecule collection"};
    }
}

auto validate_target_conformer_index(const core::MoleculeCollection& collection,
                                     const ChargeTarget& target, const std::size_t charge_set_index,
                                     const std::size_t assignment_index) -> void {
    if (!target.conformer_index.has_value()) {
        return;
    }

    const auto& molecule = collection[target.molecule_index];

    if (*target.conformer_index >= molecule.conformer_count()) {
        throw std::invalid_argument{make_context(charge_set_index, assignment_index) +
                                    "target conformer index is outside the molecule"};
    }
}

auto validate_charge_count(const core::MoleculeCollection& collection,
                           const ChargeAssignment& assignment, const std::size_t charge_set_index,
                           const std::size_t assignment_index) -> void {
    const auto& molecule = collection[assignment.target.molecule_index];

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

auto validate_assignment_order(const core::MoleculeCollection& molecules,
                               const ChargeSet& charge_set, const std::size_t charge_set_index)
    -> void {
    auto assignment_index = std::size_t{0};
    for (std::size_t molecule_index = 0; molecule_index < molecules.size(); ++molecule_index) {
        if (assignment_index == charge_set.size()) {
            throw std::invalid_argument{make_context(charge_set_index, assignment_index) +
                                        "missing target molecule"};
        }

        const auto& first_target = charge_set.assignment(assignment_index).target;
        if (first_target.molecule_index != molecule_index) {
            throw std::invalid_argument{make_context(charge_set_index, assignment_index) +
                                        "target molecule order does not match the collection"};
        }

        if (!first_target.conformer_index.has_value()) {
            ++assignment_index;
            continue;
        }

        for (std::size_t conformer_index = 0;
             conformer_index < molecules[molecule_index].conformer_count(); ++conformer_index) {
            if (assignment_index == charge_set.size()) {
                throw std::invalid_argument{make_context(charge_set_index, assignment_index) +
                                            "missing target conformer"};
            }
            const auto& target = charge_set.assignment(assignment_index).target;
            if (target.molecule_index != molecule_index ||
                target.conformer_index != conformer_index) {
                throw std::invalid_argument{make_context(charge_set_index, assignment_index) +
                                            "target conformer order does not match the molecule"};
            }
            ++assignment_index;
        }
    }

    if (assignment_index != charge_set.size()) {
        throw std::invalid_argument{make_context(charge_set_index, assignment_index) +
                                    "unexpected extra target"};
    }
}

} // namespace

auto validate_charge_collection(const core::MoleculeCollection& collection,
                                const ChargeCollection& charges) -> void {
    const auto charge_sets = charges.charge_sets();

    for (std::size_t charge_set_index = 0; charge_set_index < charge_sets.size();
         ++charge_set_index) {
        const auto assignments = charge_sets[charge_set_index].assignments();
        validate_assignment_order(collection, charge_sets[charge_set_index], charge_set_index);

        for (std::size_t assignment_index = 0; assignment_index < assignments.size();
             ++assignment_index) {
            validate_assignment(collection, assignments[assignment_index], charge_set_index,
                                assignment_index);
        }
    }
}

} // namespace chargefw::charges
