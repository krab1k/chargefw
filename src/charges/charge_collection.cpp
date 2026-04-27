#include <chargefw/charges/charge_collection.h>

#include <stdexcept>
#include <utility>

namespace chargefw::charges {
namespace {

auto validate_method_id(std::string_view method_id) -> void {
    if (method_id.empty()) {
        throw std::invalid_argument{"method id must not be empty"};
    }
}

} // namespace

ChargeSet::ChargeSet(std::string method_id, std::vector<ChargeAssignment> assignments,
                     std::optional<std::string> parameter_set_id)
    : method_id_{std::move(method_id)}, parameter_set_id_{std::move(parameter_set_id)},
      assignments_{std::move(assignments)} {
    validate_method_id(method_id_);
}

auto ChargeSet::method_id() const noexcept -> std::string_view {
    return method_id_;
}

auto ChargeSet::parameter_set_id() const noexcept -> std::optional<std::string_view> {
    if (!parameter_set_id_.has_value()) {
        return std::nullopt;
    }

    return std::string_view{*parameter_set_id_};
}

auto ChargeSet::assignments() const noexcept -> std::span<const ChargeAssignment> {
    return assignments_;
}

auto ChargeSet::assignment(const std::size_t index) const -> const ChargeAssignment& {
    return assignments_.at(index);
}

auto ChargeSet::assignment_count() const noexcept -> std::size_t {
    return assignments_.size();
}

auto ChargeSet::empty() const noexcept -> bool {
    return assignments_.empty();
}

ChargeCollection::ChargeCollection(std::vector<ChargeSet> charge_sets)
    : charge_sets_{std::move(charge_sets)} {}

auto ChargeCollection::charge_sets() const noexcept -> std::span<const ChargeSet> {
    return charge_sets_;
}

auto ChargeCollection::charge_set_count() const noexcept -> std::size_t {
    return charge_sets_.size();
}

auto ChargeCollection::empty() const noexcept -> bool {
    return charge_sets_.empty();
}

auto ChargeCollection::operator[](std::size_t index) const noexcept -> const ChargeSet& {
    return charge_sets_[index];
}

auto ChargeCollection::at(std::size_t index) const -> const ChargeSet& {
    return charge_sets_.at(index);
}

} // namespace chargefw::charges