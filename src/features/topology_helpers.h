#pragma once

#include <optional>
#include <vector>

namespace chargefw::features {

auto breadth_first_distances(const std::vector<std::vector<std::size_t>>& adjacency,
                             std::size_t source, std::optional<std::size_t> max_distance)
    -> std::vector<int>;

}
