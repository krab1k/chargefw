#include <chargefw/features/spatial_fragment.h>

#include <chargefw/features/conformer_features.h>

#include <nanoflann.hpp>

#include <algorithm>
#include <array>
#include <cmath>
#include <limits>
#include <memory>
#include <stdexcept>
#include <string>
#include <utility>
#include <vector>

namespace chargefw::features {
namespace {

inline constexpr auto no_local_index = std::numeric_limits<std::size_t>::max();

auto validate_radius(const double radius) -> void {
    if (!std::isfinite(radius) || radius <= 0.0) {
        throw std::invalid_argument{"fragment radius must be finite and strictly positive"};
    }
}

class ConformerPointCloud {
  public:
    explicit ConformerPointCloud(const ConformerFeatures& geometry) : geometry_{&geometry} {
        for (std::size_t atom_index = 0; atom_index < geometry.molecule().atom_count();
             ++atom_index) {
            const auto& position = geometry.position(atom_index);
            if (!std::isfinite(position.x) || !std::isfinite(position.y) ||
                !std::isfinite(position.z)) {
                throw std::invalid_argument{"spatial index requires finite coordinates"};
            }
        }
    }

    [[nodiscard]] auto kdtree_get_point_count() const -> std::size_t {
        return geometry_->molecule().atom_count();
    }

    [[nodiscard]] auto kdtree_get_pt(const std::size_t point_index,
                                     const std::size_t dimension) const -> double {
        const auto& position = geometry_->position(point_index);
        switch (dimension) {
        case 0:
            return position.x;
        case 1:
            return position.y;
        default:
            return position.z;
        }
    }

    template <class BoundingBox> auto kdtree_get_bbox(BoundingBox&) const -> bool {
        return false;
    }

  private:
    const ConformerFeatures* geometry_;
};

using KdTree =
    nanoflann::KDTreeSingleIndexAdaptor<nanoflann::L2_Simple_Adaptor<double, ConformerPointCloud>,
                                        ConformerPointCloud, 3, std::size_t>;

class SourceIndexRadiusResultSet {
  public:
    using DistanceType = double;
    using IndexType = std::size_t;

    SourceIndexRadiusResultSet(const double radius_squared,
                               std::vector<std::size_t>& neighbor_indices)
        : radius_squared_{radius_squared}, neighbor_indices_{&neighbor_indices} {
        neighbor_indices_->clear();
    }

    [[nodiscard]] auto size() const noexcept -> std::size_t {
        return neighbor_indices_->size();
    }
    [[nodiscard]] auto full() const noexcept -> bool {
        return true;
    }
    [[nodiscard]] auto worstDist() const noexcept -> double {
        return radius_squared_;
    }

    auto addPoint(const double squared_distance, const std::size_t point_index) -> bool {
        if (squared_distance < radius_squared_) {
            neighbor_indices_->push_back(point_index);
        }
        return true;
    }

    auto sort() -> void {
        std::ranges::sort(*neighbor_indices_);
    }

  private:
    double radius_squared_;
    std::vector<std::size_t>* neighbor_indices_;
};

} // namespace

class SpatialFragmentBuilder::SpatialIndex {
  public:
    explicit SpatialIndex(const ConformerFeatures& geometry)
        : points_{geometry}, tree_{3, points_} {}

    [[nodiscard]] auto neighbor_indices_within(const core::Position& center,
                                               const double radius) const
        -> std::vector<std::size_t> {
        const std::array<double, 3> query{center.x, center.y, center.z};
        const auto radius_squared = radius * radius;
        const auto inclusive_radius_squared =
            std::nextafter(radius_squared, std::numeric_limits<double>::infinity());

        std::vector<std::size_t> neighbors;
        SourceIndexRadiusResultSet result_set{inclusive_radius_squared, neighbors};
        static_cast<void>(tree_.radiusSearchCustomCallback(
            query.data(), result_set, nanoflann::SearchParameters{0.0F, false}));
        return neighbors;
    }

  private:
    ConformerPointCloud points_;
    KdTree tree_;
};

SpatialFragment::SpatialFragment(core::Molecule molecule,
                                 std::vector<std::size_t> local_to_source_atom_indices,
                                 std::vector<std::size_t> local_to_source_bond_indices,
                                 const std::size_t center_local_atom_index,
                                 const std::size_t source_conformer_index)
    : molecule_{std::move(molecule)},
      local_to_source_atom_indices_{std::move(local_to_source_atom_indices)},
      local_to_source_bond_indices_{std::move(local_to_source_bond_indices)},
      center_local_atom_index_{center_local_atom_index},
      source_conformer_index_{source_conformer_index} {}

auto SpatialFragment::molecule() const noexcept -> const core::Molecule& {
    return molecule_;
}

auto SpatialFragment::source_conformer_index() const noexcept -> std::size_t {
    return source_conformer_index_;
}

auto SpatialFragment::center_local_atom_index() const noexcept -> std::size_t {
    return center_local_atom_index_;
}

auto SpatialFragment::local_to_source_atom_indices() const noexcept
    -> std::span<const std::size_t> {
    return local_to_source_atom_indices_;
}

auto SpatialFragment::local_to_source_bond_indices() const noexcept
    -> std::span<const std::size_t> {
    return local_to_source_bond_indices_;
}

SpatialFragmentBuilder::SpatialFragmentBuilder(const PreparedMolecule& source,
                                               const ConformerFeatures& geometry)
    : source_{&source}, geometry_{&geometry} {
    if (std::addressof(geometry.molecule()) != std::addressof(source.molecule())) {
        throw std::invalid_argument{"fragment geometry must refer to the source molecule"};
    }
    spatial_index_ = std::make_unique<SpatialIndex>(geometry);
}

SpatialFragmentBuilder::~SpatialFragmentBuilder() = default;

auto SpatialFragmentBuilder::atom_indices_within(const std::size_t center_atom_index,
                                                 const double radius) const
    -> std::vector<std::size_t> {
    validate_radius(radius);

    auto atom_indices =
        spatial_index_->neighbor_indices_within(geometry_->position(center_atom_index), radius);
    std::ranges::sort(atom_indices);
    return atom_indices;
}

auto SpatialFragmentBuilder::build(const std::size_t center_atom_index, const double radius) const
    -> SpatialFragment {
    const auto& molecule = source_->molecule();
    const auto& topology = source_->topology();
    const auto& geometry = *geometry_;

    auto selected_source_atom_indices = atom_indices_within(center_atom_index, radius);
    // Local atoms follow sorted source order, so a source index maps to its local position by
    // binary search. This keeps fragment construction proportional to the fragment size.
    const auto local_index_of = [&](const std::size_t source_index) -> std::size_t {
        const auto found = std::ranges::lower_bound(selected_source_atom_indices, source_index);
        return found != selected_source_atom_indices.end() && *found == source_index
                   ? static_cast<std::size_t>(found - selected_source_atom_indices.begin())
                   : no_local_index;
    };

    std::vector<core::Atom> atoms;
    std::vector<core::Position> positions;
    atoms.reserve(selected_source_atom_indices.size());
    positions.reserve(selected_source_atom_indices.size());

    auto center_local_atom_index = no_local_index;
    for (std::size_t local_index = 0; local_index < selected_source_atom_indices.size();
         ++local_index) {
        const auto source_index = selected_source_atom_indices[local_index];
        const auto& source_atom = molecule.atom(source_index);
        atoms.emplace_back(source_atom.atomic_number(), source_atom.formal_charge(),
                           std::string{source_atom.name()});
        positions.push_back(geometry.position(source_index));

        if (source_index == center_atom_index) {
            center_local_atom_index = local_index;
        }
    }

    // Each bond with both atoms selected appears in the incident lists of both atoms; keep it once.
    // Sorting restores source bond order.
    std::vector<std::size_t> local_to_source_bond_indices;
    for (const auto source_index : selected_source_atom_indices) {
        for (const auto source_bond_index : topology.incident_bond_indices(source_index)) {
            const auto& source_bond = molecule.bond(source_bond_index);
            const auto other_index = source_bond.first_atom_index() == source_index
                                         ? source_bond.second_atom_index()
                                         : source_bond.first_atom_index();
            if (source_index < other_index && local_index_of(other_index) != no_local_index) {
                local_to_source_bond_indices.push_back(source_bond_index);
            }
        }
    }
    std::ranges::sort(local_to_source_bond_indices);

    std::vector<core::Bond> bonds;
    bonds.reserve(local_to_source_bond_indices.size());
    for (const auto source_bond_index : local_to_source_bond_indices) {
        const auto& source_bond = molecule.bond(source_bond_index);
        bonds.emplace_back(local_index_of(source_bond.first_atom_index()),
                           local_index_of(source_bond.second_atom_index()), source_bond.order());
    }

    std::vector<core::Conformer> conformers;
    conformers.emplace_back(std::move(positions),
                            std::string{molecule.conformer(geometry.conformer_index()).name()});
    return SpatialFragment{core::Molecule{std::move(atoms), std::move(bonds), std::move(conformers),
                                          std::string{molecule.name()}},
                           std::move(selected_source_atom_indices),
                           std::move(local_to_source_bond_indices), center_local_atom_index,
                           geometry.conformer_index()};
}

} // namespace chargefw::features
