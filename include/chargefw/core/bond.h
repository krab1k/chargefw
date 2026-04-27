#pragma once

#include <cstddef>

namespace chargefw::core {
enum class BondOrder { UNKNOWN, SINGLE, DOUBLE, TRIPLE, AROMATIC };

class Bond {
  public:
    explicit Bond(std::size_t first_atom_index, std::size_t second_atom_index,
                  BondOrder order = BondOrder::SINGLE);

    [[nodiscard]] auto first_atom_index() const noexcept -> std::size_t;
    [[nodiscard]] auto second_atom_index() const noexcept -> std::size_t;
    [[nodiscard]] auto order() const noexcept -> BondOrder;

  private:
    std::size_t first_atom_index_;
    std::size_t second_atom_index_;
    BondOrder order_;
};

} // namespace chargefw::core