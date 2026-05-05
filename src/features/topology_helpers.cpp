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

} // namespace chargefw::features
