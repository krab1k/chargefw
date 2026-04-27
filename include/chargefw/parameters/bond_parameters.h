#pragma once

#include <chargefw/parameters/atom_parameters.h>
#include <chargefw/parameters/common_parameters.h>

#include <cstddef>
#include <span>
#include <string>
#include <string_view>
#include <vector>

namespace chargefw::parameters {

struct BondParameterKey {
    AtomParameterKey first_atom;
    AtomParameterKey second_atom;
    std::string bond_class;
    std::string bond_type;
};

struct BondParameterEntry {
    BondParameterKey key;
    std::vector<NamedParameter> parameters;
};

class BondParameters {
public:
    BondParameters() = default;

    explicit BondParameters(std::vector<BondParameterEntry> entries);

    [[nodiscard]] auto entries() const noexcept -> std::span<const BondParameterEntry>;

    [[nodiscard]] auto size() const noexcept -> std::size_t;
    [[nodiscard]] auto empty() const noexcept -> bool;

    [[nodiscard]] auto entry(std::size_t index) const -> const BondParameterEntry&;

    [[nodiscard]] auto operator[](std::size_t index) const noexcept -> const BondParameterEntry&;
    [[nodiscard]] auto at(std::size_t index) const -> const BondParameterEntry&;

    [[nodiscard]] auto contains(std::size_t entry_index, std::string_view name) const -> bool;

    [[nodiscard]] auto parameter(std::size_t entry_index, std::string_view name) const -> double;

private:
    std::vector<BondParameterEntry> entries_;
};

} // namespace chargefw::parameters