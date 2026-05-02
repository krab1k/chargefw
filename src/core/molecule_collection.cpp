#include <chargefw/core/molecule_collection.h>

#include <utility>

namespace chargefw::core {

MoleculeCollection::MoleculeCollection(std::vector<Molecule> molecules, std::string name)
    : molecules_{std::move(molecules)}, name_{std::move(name)} {}

auto MoleculeCollection::name() const noexcept -> std::string_view {
    return name_;
}

auto MoleculeCollection::molecules() const noexcept -> std::span<const Molecule> {
    return molecules_;
}

auto MoleculeCollection::molecule_count() const noexcept -> std::size_t {
    return molecules_.size();
}

auto MoleculeCollection::empty() const noexcept -> bool {
    return molecules_.empty();
}

[[nodiscard]] auto MoleculeCollection::operator[](std::size_t index) const noexcept
    -> const Molecule& {
    return molecules_[index];
}

[[nodiscard]] auto MoleculeCollection::at(std::size_t index) const -> const Molecule& {
    return molecules_.at(index);
}

} // namespace chargefw::core