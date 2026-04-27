#include <chargefw/methods/method_options.h>

#include <utility>

namespace chargefw::methods {

MethodOptions::MethodOptions(std::unordered_map<std::string, MethodOptionValue> values)
    : values_{std::move(values)} {}

auto MethodOptions::contains(const std::string_view id) const -> bool {
    return values_.contains(std::string{id});
}

auto MethodOptions::set(std::string id, MethodOptionValue value) -> void {
    values_.insert_or_assign(std::move(id), std::move(value));
}

auto make_default_options(std::span<const MethodOptionSpec> schema) -> MethodOptions {
    MethodOptions options;

    for (const auto& option : schema) {
        options.set(std::string{option.id}, option.default_value);
    }

    return options;
}

auto validate_method_options(std::span<const MethodOptionSpec> schema, const MethodOptions& options)
    -> void {
    for (const auto& option : schema) {
        if (!options.contains(option.id)) {
            throw std::invalid_argument{"missing method option '" + std::string{option.id} + "'"};
        }
    }
}

} // namespace chargefw::methods