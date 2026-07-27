#pragma once

#include <optional>
#include <vector>

namespace chargefw::features {

auto breadth_first_distances(const std::vector<std::vector<std::size_t>>& adjacency,
                             std::size_t source, std::optional<std::size_t> max_distance)
    -> std::vector<int>;

auto all_pairs_bond_distances(const std::vector<std::vector<std::size_t>>& adjacency)
    -> std::vector<std::vector<int>>;

auto is_connected(const std::vector<std::vector<std::size_t>>& adjacency) -> bool;

} // namespace chargefw::features
