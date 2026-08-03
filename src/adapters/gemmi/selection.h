#pragma once

#include <chargefw/adapters/gemmi/input_options.h>

#include <gemmi/model.hpp>

#include <cstddef>
#include <optional>
#include <string>
#include <string_view>
#include <unordered_map>
#include <utility>
#include <vector>

namespace chargefw::adapters::gemmi::selection {

struct SelectedResidue {
    const ::gemmi::Residue* residue;
    std::vector<std::pair<std::string_view, std::size_t>> atom_indices;

    [[nodiscard]] auto find_atom(std::string_view name) const -> std::optional<std::size_t>;
};

// Borrows a Gemmi model and its atoms. The model must outlive this view and remain unmodified.
class SelectedModel {
  public:
    explicit SelectedModel(const ::gemmi::Model& model, RecordSelection selection);

    [[nodiscard]] auto atoms() const -> const std::vector<const ::gemmi::Atom*>&;
    [[nodiscard]] auto residues() const -> const std::vector<SelectedResidue>&;
    [[nodiscard]] auto atom_index(const ::gemmi::Atom* atom) const -> std::optional<std::size_t>;
    [[nodiscard]] auto atom_index_by_serial(int serial) const -> std::optional<std::size_t>;

  private:
    std::vector<const ::gemmi::Atom*> atoms_;
    std::vector<SelectedResidue> residues_;
    std::unordered_map<const ::gemmi::Atom*, std::size_t> atom_indices_;
    std::unordered_map<int, std::size_t> serial_indices_;
};

} // namespace chargefw::adapters::gemmi::selection
