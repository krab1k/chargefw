#include "support/test_assertions.h"

#include <chargefw/methods/method_options.h>

#include <array>
#include <cassert>
#include <stdexcept>
#include <string>
#include <string_view>

namespace methods = chargefw::methods;

auto main() -> int {
    const std::array schema{
        methods::MethodOptionSpec{.id = "enabled",
                                  .description = "Enable the method",
                                  .type = methods::MethodOptionType::boolean,
                                  .default_value = true,
                                  .choices = {}},
        methods::MethodOptionSpec{.id = "limit",
                                  .description = "Iteration limit",
                                  .type = methods::MethodOptionType::integer,
                                  .default_value = 25,
                                  .choices = {}},
        methods::MethodOptionSpec{.id = "radius",
                                  .description = "Search radius",
                                  .type = methods::MethodOptionType::floating_point,
                                  .default_value = 1.5,
                                  .choices = {},
                                  .minimum = 0.0,
                                  .minimum_inclusive = false},
        methods::MethodOptionSpec{.id = "label",
                                  .description = "Label",
                                  .type = methods::MethodOptionType::string,
                                  .default_value = std::string{"default"},
                                  .choices = {}}};

    auto options = methods::make_default_options(schema);

    assert(options.contains("enabled"));
    assert(options.contains("limit"));
    assert(options.contains("radius"));
    assert(options.contains("label"));
    assert(options.get<bool>("enabled"));
    assert(options.get<int>("limit") == 25);
    assert(options.get<double>("radius") == 1.5);
    assert(options.get<std::string>("label") == "default");

    options.set("limit", 50);
    assert(options.get<int>("limit") == 50);

    methods::validate_method_options(schema, options);

    assert(chargefw::test::throws_invalid_argument(
        [&schema] -> void { methods::validate_method_options(schema, methods::MethodOptions{}); }));

    assert(chargefw::test::throws_invalid_argument([&schema] -> void {
        auto invalid_options = methods::make_default_options(schema);
        invalid_options.set("limit", 50.0);
        methods::validate_method_options(schema, invalid_options);
    }));

    const std::array mismatched_bound_schema{
        methods::MethodOptionSpec{.id = "limit",
                                  .description = "Iteration limit",
                                  .type = methods::MethodOptionType::integer,
                                  .default_value = 1,
                                  .choices = {},
                                  .minimum = 0.0}};
    assert(chargefw::test::throws_invalid_argument([&mismatched_bound_schema] -> void {
        methods::validate_method_option_schema(mismatched_bound_schema);
    }));

    assert(chargefw::test::throws_invalid_argument([&schema] -> void {
        auto invalid_options = methods::make_default_options(schema);
        invalid_options.set("radius", 0.0);
        methods::validate_method_options(schema, invalid_options);
    }));

    assert(chargefw::test::throws_invalid_argument([&schema] -> void {
        auto invalid_options = methods::make_default_options(schema);
        invalid_options.set("extra", true);
        methods::validate_method_options(schema, invalid_options);
    }));

    bool rejected_unknown_option_lookup = false;

    try {
        [[maybe_unused]] const auto& invalid = options.get<int>("missing");
    } catch (const std::out_of_range&) {
        rejected_unknown_option_lookup = true;
    }

    assert(rejected_unknown_option_lookup);

    bool rejected_wrong_option_type = false;

    try {
        [[maybe_unused]] const auto& invalid = options.get<double>("limit");
    } catch (const std::invalid_argument&) {
        rejected_wrong_option_type = true;
    }

    assert(rejected_wrong_option_type);

    const std::array choice_schema{
        methods::MethodOptionSpec{.id = "mode",
                                  .description = "Calculation mode",
                                  .type = methods::MethodOptionType::string,
                                  .default_value = std::string{"fast"},
                                  .choices = {std::string{"fast"}, std::string{"accurate"}}},
        methods::MethodOptionSpec{.id = "limit",
                                  .description = "Iteration limit",
                                  .type = methods::MethodOptionType::integer,
                                  .default_value = 25,
                                  .choices = {10, 25, 50}}};

    auto choice_options = methods::make_default_options(choice_schema);
    methods::validate_method_options(choice_schema, choice_options);

    assert(chargefw::test::throws_invalid_argument([&choice_schema] -> void {
        auto invalid_options = methods::make_default_options(choice_schema);
        invalid_options.set("mode", std::string{"slow"});
        methods::validate_method_options(choice_schema, invalid_options);
    }));

    assert(chargefw::test::throws_invalid_argument([&choice_schema] -> void {
        auto invalid_options = methods::make_default_options(choice_schema);
        invalid_options.set("limit", 100);
        methods::validate_method_options(choice_schema, invalid_options);
    }));

    const std::array duplicate_id_schema{
        methods::MethodOptionSpec{.id = "limit",
                                  .description = "First limit",
                                  .type = methods::MethodOptionType::integer,
                                  .default_value = 25,
                                  .choices = {}},
        methods::MethodOptionSpec{.id = "limit",
                                  .description = "Second limit",
                                  .type = methods::MethodOptionType::integer,
                                  .default_value = 50,
                                  .choices = {}}};

    assert(chargefw::test::throws_invalid_argument([&duplicate_id_schema] -> void {
        methods::validate_method_option_schema(duplicate_id_schema);
    }));

    const std::array wrong_default_type_schema{
        methods::MethodOptionSpec{.id = "limit",
                                  .description = "Iteration limit",
                                  .type = methods::MethodOptionType::integer,
                                  .default_value = 25.0,
                                  .choices = {}}};

    assert(chargefw::test::throws_invalid_argument([&wrong_default_type_schema] -> void {
        methods::validate_method_option_schema(wrong_default_type_schema);
    }));

    const std::array wrong_choice_type_schema{
        methods::MethodOptionSpec{.id = "limit",
                                  .description = "Iteration limit",
                                  .type = methods::MethodOptionType::integer,
                                  .default_value = 25,
                                  .choices = {25, 50.0}}};

    assert(chargefw::test::throws_invalid_argument([&wrong_choice_type_schema] -> void {
        methods::validate_method_option_schema(wrong_choice_type_schema);
    }));

    const std::array default_not_allowed_schema{
        methods::MethodOptionSpec{.id = "limit",
                                  .description = "Iteration limit",
                                  .type = methods::MethodOptionType::integer,
                                  .default_value = 25,
                                  .choices = {50, 100}}};

    assert(chargefw::test::throws_invalid_argument([&default_not_allowed_schema] -> void {
        methods::validate_method_option_schema(default_not_allowed_schema);
    }));

    return 0;
}
