#include <algorithm>
#include <optional>
#include <queue>
#include <utility>
#include <vector>

namespace chargefw::features {

auto breadth_first_distances(const std::vector<std::vector<std::size_t>>& adjacency,
                             const std::size_t source,
                             const std::optional<std::size_t> max_distance) -> std::vector<int> {
    auto distances = std::vector(adjacency.size(), -1);
    auto queue = std::queue<std::size_t>{};

    distances[source] = 0;
    queue.push(source);

    while (!queue.empty()) {
        const auto current = queue.front();
        queue.pop();

        const auto current_distance = distances[current];

        if (max_distance.has_value() && std::cmp_greater_equal(current_distance, *max_distance)) {
            continue;
        }

        for (const auto neighbor : adjacency[current]) {
            if (distances[neighbor] != -1) {
                continue;
            }

            distances[neighbor] = current_distance + 1;
            queue.push(neighbor);
        }
    }

    return distances;
}

auto all_pairs_bond_distances(const std::vector<std::vector<std::size_t>>& adjacency)
    -> std::vector<std::vector<int>> {
    auto distances = std::vector<std::vector<int>>{};
    distances.reserve(adjacency.size());

    for (std::size_t atom_index = 0; atom_index < adjacency.size(); ++atom_index) {
        distances.push_back(breadth_first_distances(adjacency, atom_index, std::nullopt));
    }

    return distances;
}

auto is_connected(const std::vector<std::vector<std::size_t>>& adjacency) -> bool {
    if (adjacency.empty()) {
        return true;
    }

    const auto distances = breadth_first_distances(adjacency, 0, std::nullopt);

    return std::ranges::all_of(distances, [](const int distance) -> bool { return distance >= 0; });
}

} // namespace chargefw::features
