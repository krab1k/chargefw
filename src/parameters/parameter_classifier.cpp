#include <chargefw/parameters/parameter_classifier.h>

#include <algorithm>
#include <stdexcept>
#include <string>
#include <vector>

namespace chargefw::parameters {
namespace {

auto bond_order_type(const core::BondOrder order) -> std::string
{
    switch (order) {
    case core::BondOrder::UNKNOWN:
        return "unknown";

    case core::BondOrder::SINGLE:
        return "1";

    case core::BondOrder::DOUBLE:
        return "2";

    case core::BondOrder::TRIPLE:
        return "3";

    case core::BondOrder::AROMATIC:
        return "aromatic";
    }

    return "unknown";
}

auto bond_order_rank(const core::BondOrder order) -> int
{
    switch (order) {
    case core::BondOrder::UNKNOWN:
        return 0;

    case core::BondOrder::SINGLE:
        return 1;

    case core::BondOrder::DOUBLE:
        return 2;

    case core::BondOrder::TRIPLE:
        return 3;

    case core::BondOrder::AROMATIC:
        return 2;
    }

    return 0;
}

auto highest_bond_order_type(const core::Molecule& molecule,
                             const features::TopologyFeatures& topology,
                             const std::size_t atom_index) -> std::string
{
    auto highest = 0;

    for (const auto bond_index : topology.incident_bond_indices(atom_index)) {
        highest = std::max(highest, bond_order_rank(molecule.bond(bond_index).order()));
    }

    return std::to_string(highest);
}

auto atom_type_for(const core::Molecule& molecule,
                   const features::TopologyFeatures& topology,
                   const std::size_t atom_index,
                   const AtomParameterClassificationKind classification) -> std::string
{
    switch (classification) {
    case AtomParameterClassificationKind::PLAIN:
        return "*";

    case AtomParameterClassificationKind::HIGHEST_BOND_ORDER:
        return highest_bond_order_type(molecule, topology, atom_index);

    case AtomParameterClassificationKind::BONDED_ELEMENTS:
        throw std::logic_error{
            "bonded-elements atom parameter classification is not implemented yet"
        };
    }

    throw std::logic_error{"unknown atom parameter classification"};
}

auto bond_type_for(const core::Molecule& molecule,
                   const std::size_t bond_index,
                   const BondParameterClassificationKind classification) -> std::string
{
    switch (classification) {
    case BondParameterClassificationKind::PLAIN:
        return "*";

    case BondParameterClassificationKind::BOND_ORDER:
        return bond_order_type(molecule.bond(bond_index).order());
    }

    throw std::logic_error{"unknown bond parameter classification"};
}

auto matches_atom_key(const core::Molecule& molecule,
                      const features::TopologyFeatures& topology,
                      const std::size_t atom_index,
                      const AtomParameterKey& key) -> bool
{
    const auto& atom = molecule.atom(atom_index);

    if (key.atomic_number != 0 && key.atomic_number != atom.atomic_number()) {
        return false;
    }

    return key.type == atom_type_for(molecule, topology, atom_index, key.classification);
}

auto find_atom_parameter_entry(const core::Molecule& molecule,
                               const features::TopologyFeatures& topology,
                               const std::size_t atom_index,
                               const AtomParameters& parameters) -> std::size_t
{
    for (std::size_t entry_index = 0; entry_index < parameters.size(); ++entry_index) {
        const auto& [key, named_parameters] = parameters[entry_index];

        if (matches_atom_key(molecule, topology, atom_index, key)) {
            return entry_index;
        }
    }

    throw std::invalid_argument{
        "no atom parameter entry for atom index " + std::to_string(atom_index)
    };
}

auto classify_atoms(const core::Molecule& molecule,
                    const features::TopologyFeatures& topology,
                    const AtomParameters& parameters) -> AtomParameterClassification
{
    if (parameters.empty() && molecule.atom_count() != 0) {
        throw std::invalid_argument{"atom parameters are required but empty"};
    }

    std::vector<std::size_t> indices;
    indices.reserve(molecule.atom_count());

    for (std::size_t atom_index = 0; atom_index < molecule.atom_count(); ++atom_index) {
        indices.push_back(find_atom_parameter_entry(molecule, topology, atom_index, parameters));
    }

    return AtomParameterClassification{std::move(indices)};
}

auto matches_bond_key(const core::Molecule& molecule,
                      const features::TopologyFeatures& topology,
                      const std::size_t bond_index,
                      const BondParameterKey& key) -> bool
{
    const auto& bond = molecule.bond(bond_index);

    const auto first_atom_index = bond.first_atom_index();
    const auto second_atom_index = bond.second_atom_index();

    const auto forward =
        matches_atom_key(molecule, topology, first_atom_index, key.first_atom) &&
        matches_atom_key(molecule, topology, second_atom_index, key.second_atom);

    const auto reverse =
        matches_atom_key(molecule, topology, first_atom_index, key.second_atom) &&
        matches_atom_key(molecule, topology, second_atom_index, key.first_atom);

    if (!forward && !reverse) {
        return false;
    }

    return key.bond.type == bond_type_for(molecule, bond_index, key.bond.classification);
}

auto find_bond_parameter_entry(const core::Molecule& molecule,
                               const features::TopologyFeatures& topology,
                               const std::size_t bond_index,
                               const BondParameters& parameters) -> std::size_t
{
    for (std::size_t entry_index = 0; entry_index < parameters.size(); ++entry_index) {
        const auto& [key, named_parameters] = parameters[entry_index];

        if (matches_bond_key(molecule, topology, bond_index, key)) {
            return entry_index;
        }
    }

    throw std::invalid_argument{
        "no bond parameter entry for bond index " + std::to_string(bond_index)
    };
}

auto classify_bonds(const core::Molecule& molecule,
                    const features::TopologyFeatures& topology,
                    const BondParameters& parameters) -> BondParameterClassification
{
    if (parameters.empty()) {
        return BondParameterClassification{};
    }

    std::vector<std::size_t> indices;
    indices.reserve(molecule.bond_count());

    for (std::size_t bond_index = 0; bond_index < molecule.bond_count(); ++bond_index) {
        indices.push_back(find_bond_parameter_entry(molecule, topology, bond_index, parameters));
    }

    return BondParameterClassification{std::move(indices)};
}

} // namespace

auto classify_parameters(const core::Molecule& molecule,
                         const features::TopologyFeatures& topology,
                         const ParameterSet& parameters) -> ParameterClassification
{
    auto atom_classification = classify_atoms(molecule, topology, parameters.atom());
    auto bond_classification = classify_bonds(molecule, topology, parameters.bond());

    auto classification = ParameterClassification{
        std::move(atom_classification),
        std::move(bond_classification)
    };

    validate_parameter_classification(molecule, parameters, classification);

    return classification;
}

} // namespace chargefw::parameters