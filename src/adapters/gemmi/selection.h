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

[[nodiscard]] auto include_residue(const ::gemmi::Residue& residue, RecordSelection selection)
    -> bool;

[[nodiscard]] auto is_first_named_atom(const ::gemmi::Residue& residue, std::size_t atom_index)
    -> bool;

[[nodiscard]] auto select_altloc(const ::gemmi::Residue& residue, std::size_t first_atom_index)
    -> const ::gemmi::Atom&;

struct SelectedResidue {
    const ::gemmi::Residue* residue;
    std::vector<std::pair<std::string, std::size_t>> atom_indices;
};

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

[[nodiscard]] auto atom_index(const SelectedResidue& residue, std::string_view name)
    -> std::optional<std::size_t>;

} // namespace chargefw::adapters::gemmi::selection
