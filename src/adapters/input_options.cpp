#include <chargefw/adapters/conformer_selection.h>
#include <chargefw/adapters/gemmi/input_options.h>

#include <stdexcept>
#include <string>

namespace chargefw::adapters {

auto conformer_selection_from_string(const std::string_view value) -> ConformerSelection {
    if (value == "first") {
        return ConformerSelection::first;
    }
    if (value == "all") {
        return ConformerSelection::all;
    }
    throw std::invalid_argument{"unknown conformer selection: " + std::string{value}};
}

auto to_string(const ConformerSelection value) -> std::string_view {
    switch (value) {
    case ConformerSelection::first:
        return "first";
    case ConformerSelection::all:
        return "all";
    }
    throw std::invalid_argument{"unknown conformer selection"};
}

} // namespace chargefw::adapters

namespace chargefw::adapters::gemmi {

auto record_selection_from_string(const std::string_view value) -> RecordSelection {
    if (value == "all") {
        return RecordSelection::all;
    }
    if (value == "polymers-and-ligands") {
        return RecordSelection::polymers_and_ligands;
    }
    if (value == "polymers") {
        return RecordSelection::polymers;
    }
    throw std::invalid_argument{"unknown record selection: " + std::string{value}};
}

auto to_string(const RecordSelection value) -> std::string_view {
    switch (value) {
    case RecordSelection::all:
        return "all";
    case RecordSelection::polymers_and_ligands:
        return "polymers-and-ligands";
    case RecordSelection::polymers:
        return "polymers";
    }
    throw std::invalid_argument{"unknown record selection"};
}

auto bond_strategy_from_string(const std::string_view value) -> BondStrategy {
    if (value == "none") {
        return BondStrategy::none;
    }
    if (value == "explicit") {
        return BondStrategy::explicit_bonds;
    }
    if (value == "templates") {
        return BondStrategy::templates;
    }
    if (value == "hybrid") {
        return BondStrategy::hybrid;
    }
    throw std::invalid_argument{"unknown bond strategy: " + std::string{value}};
}

auto to_string(const BondStrategy value) -> std::string_view {
    switch (value) {
    case BondStrategy::none:
        return "none";
    case BondStrategy::explicit_bonds:
        return "explicit";
    case BondStrategy::templates:
        return "templates";
    case BondStrategy::hybrid:
        return "hybrid";
    }
    throw std::invalid_argument{"unknown bond strategy"};
}

} // namespace chargefw::adapters::gemmi
