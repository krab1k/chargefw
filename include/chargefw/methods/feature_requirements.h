#pragma once

namespace chargefw::methods {

struct FeatureRequirements {
    bool bond_graph = false;
    bool bond_orders = false;
    bool topological_distances = false;

    bool formal_charges = false;
    bool element_properties = false;

    bool common_parameters = false;
    bool atom_parameters = false;
    bool bond_parameters = false;

    bool coordinates = false;
    bool spatial_neighbor_search = false;
};
} // namespace chargefw::methods
