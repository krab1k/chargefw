#pragma once

#include <chargefw/core/molecule.h>

#include <cstddef>
#include <span>
#include <string>
#include <string_view>
#include <vector>

namespace chargefw::core {

class MoleculeCollection {
  public:
    explicit MoleculeCollection(std::vector<Molecule> molecules, std::string name = {});

    [[nodiscard]] auto name() const noexcept -> std::string_view;

    [[nodiscard]] auto molecules() const noexcept -> std::span<const Molecule>;

    [[nodiscard]] auto molecule_count() const noexcept -> std::size_t;

    [[nodiscard]] auto empty() const noexcept -> bool;

    [[nodiscard]] auto operator[](std::size_t index) const noexcept -> const Molecule&;
    [[nodiscard]] auto at(std::size_t index) const -> const Molecule&;

  private:
    std::vector<Molecule> molecules_;
    std::string name_;
};

} // namespace chargefw::core