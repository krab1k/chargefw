#include <chargefw/parameters/classification/classification_result.h>

#include <stdexcept>
#include <utility>

namespace chargefw::parameters {

ClassificationResult::ClassificationResult(ParameterClassification classification)
    : classification_{std::move(classification)} {}

ClassificationResult::ClassificationResult(std::vector<ClassificationIssue> issues)
    : issues_{std::move(issues)} {}

ClassificationResult::operator bool() const noexcept {
    return issues_.empty();
}

auto ClassificationResult::classification() const -> const ParameterClassification& {
    if (!issues_.empty()) {
        throw std::logic_error{"classification is not available for a failed result"};
    }

    return classification_;
}

auto ClassificationResult::issues() const noexcept -> const std::vector<ClassificationIssue>& {
    return issues_;
}

} // namespace chargefw::parameters
