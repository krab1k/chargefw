#include <chargefw/methods/method_prerequisites.h>

#include <utility>

namespace chargefw::methods {

auto PrerequisiteResult::ok() const noexcept -> bool {
    return issues_.empty();
}

PrerequisiteResult::operator bool() const noexcept {
    return ok();
}

auto PrerequisiteResult::issues() const noexcept -> std::span<const PrerequisiteIssue> {
    return issues_;
}

auto PrerequisiteResult::add(PrerequisiteIssue issue) -> void {
    issues_.push_back(std::move(issue));
}

} // namespace chargefw::methods