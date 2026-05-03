#include <chargefw/features/prepared_molecule.h>
#include <chargefw/features/prepared_molecule_collection.h>

namespace chargefw::features {

PreparedMoleculeCollection::PreparedMoleculeCollection(const core::MoleculeCollection& collection) {
    molecules_.reserve(collection.size());

    for (const auto& molecule : collection.molecules()) {
        molecules_.emplace_back(molecule);
    }
}

auto PreparedMoleculeCollection::molecules() const noexcept -> std::span<const PreparedMolecule> {
    return molecules_;
}

auto PreparedMoleculeCollection::size() const noexcept -> std::size_t {
    return molecules_.size();
}

auto PreparedMoleculeCollection::empty() const noexcept -> bool {
    return molecules_.empty();
}

auto PreparedMoleculeCollection::operator[](std::size_t index) const noexcept
    -> const PreparedMolecule& {
    return molecules_[index];
}

auto PreparedMoleculeCollection::at(std::size_t index) const -> const PreparedMolecule& {
    return molecules_.at(index);
}

} // namespace chargefw::features