#include <chargefw/parameters/parameter_set.h>

#include <stdexcept>
#include <utility>

namespace chargefw::parameters {
namespace {

auto validate_metadata(const ParameterSetMetadata& metadata) -> void {
    if (metadata.id.empty()) {
        throw std::invalid_argument{"parameter set id must not be empty"};
    }

    if (metadata.method_id.empty()) {
        throw std::invalid_argument{"parameter set method id must not be empty"};
    }
}

} // namespace

ParameterSet::ParameterSet(ParameterSetMetadata metadata, CommonParameters common,
                           AtomParameters atom, BondParameters bond)
    : metadata_{std::move(metadata)}, common_{std::move(common)}, atom_{std::move(atom)},
      bond_{std::move(bond)} {
    validate_metadata(metadata_);
}

auto ParameterSet::metadata() const noexcept -> const ParameterSetMetadata& {
    return metadata_;
}

auto ParameterSet::id() const noexcept -> std::string_view {
    return metadata_.id;
}

auto ParameterSet::method_id() const noexcept -> std::string_view {
    return metadata_.method_id;
}

auto ParameterSet::name() const noexcept -> std::string_view {
    return metadata_.name;
}

auto ParameterSet::publication() const noexcept -> std::string_view {
    return metadata_.publication;
}

auto ParameterSet::notes() const noexcept -> std::string_view {
    return metadata_.notes;
}

auto ParameterSet::common() const noexcept -> const CommonParameters& {
    return common_;
}

auto ParameterSet::atom() const noexcept -> const AtomParameters& {
    return atom_;
}

auto ParameterSet::bond() const noexcept -> const BondParameters& {
    return bond_;
}

} // namespace chargefw::parameters