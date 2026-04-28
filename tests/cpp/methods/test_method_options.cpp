#include <chargefw/methods/method_options.h>

#include <cassert>
#include <stdexcept>
#include <string>
#include <string_view>

namespace methods = chargefw::methods;

namespace {

template <typename Function>
[[nodiscard]] auto rejects_invalid_argument(Function function) -> bool {
    try {
        function();
    } catch (const std::invalid_argument&) {
        return true;
    }

    return false;
}

} // namespace

auto main() -> int {
    const methods::MethodOptionSpec schema[]{{.id = "enabled",
                                              .description = "Enable the method",
                                              .type = methods::MethodOptionType::boolean,
                                              .default_value = true,
                                              .choices = {}},
                                             {.id = "limit",
                                              .description = "Iteration limit",
                                              .type = methods::MethodOptionType::integer,
                                              .default_value = 25,
                                              .choices = {}},
                                             {.id = "radius",
                                              .description = "Search radius",
                                              .type = methods::MethodOptionType::floating_point,
                                              .default_value = 1.5,
                                              .choices = {}},
                                             {.id = "label",
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

    assert(rejects_invalid_argument(
        [&schema] { methods::validate_method_options(schema, methods::MethodOptions{}); }));

    assert(rejects_invalid_argument([&schema] {
        auto invalid_options = methods::make_default_options(schema);
        invalid_options.set("limit", 50.0);
        methods::validate_method_options(schema, invalid_options);
    }));

    assert(rejects_invalid_argument([&schema] {
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

    const methods::MethodOptionSpec choice_schema[]{
        {.id = "mode",
         .description = "Calculation mode",
         .type = methods::MethodOptionType::string,
         .default_value = std::string{"fast"},
         .choices = {std::string{"fast"}, std::string{"accurate"}}},
        {.id = "limit",
         .description = "Iteration limit",
         .type = methods::MethodOptionType::integer,
         .default_value = 25,
         .choices = {10, 25, 50}}};

    auto choice_options = methods::make_default_options(choice_schema);
    methods::validate_method_options(choice_schema, choice_options);

    assert(rejects_invalid_argument([&choice_schema] {
        auto invalid_options = methods::make_default_options(choice_schema);
        invalid_options.set("mode", std::string{"slow"});
        methods::validate_method_options(choice_schema, invalid_options);
    }));

    assert(rejects_invalid_argument([&choice_schema] {
        auto invalid_options = methods::make_default_options(choice_schema);
        invalid_options.set("limit", 100);
        methods::validate_method_options(choice_schema, invalid_options);
    }));

    const methods::MethodOptionSpec duplicate_id_schema[]{
        {.id = "limit",
         .description = "First limit",
         .type = methods::MethodOptionType::integer,
         .default_value = 25,
         .choices = {}},
        {.id = "limit",
         .description = "Second limit",
         .type = methods::MethodOptionType::integer,
         .default_value = 50,
         .choices = {}}};

    assert(rejects_invalid_argument(
        [&duplicate_id_schema] { methods::validate_method_option_schema(duplicate_id_schema); }));

    const methods::MethodOptionSpec wrong_default_type_schema[]{
        {.id = "limit",
         .description = "Iteration limit",
         .type = methods::MethodOptionType::integer,
         .default_value = 25.0,
         .choices = {}}};

    assert(rejects_invalid_argument([&wrong_default_type_schema] {
        methods::validate_method_option_schema(wrong_default_type_schema);
    }));

    const methods::MethodOptionSpec wrong_choice_type_schema[]{
        {.id = "limit",
         .description = "Iteration limit",
         .type = methods::MethodOptionType::integer,
         .default_value = 25,
         .choices = {25, 50.0}}};

    assert(rejects_invalid_argument([&wrong_choice_type_schema] {
        methods::validate_method_option_schema(wrong_choice_type_schema);
    }));

    const methods::MethodOptionSpec default_not_allowed_schema[]{
        {.id = "limit",
         .description = "Iteration limit",
         .type = methods::MethodOptionType::integer,
         .default_value = 25,
         .choices = {50, 100}}};

    assert(rejects_invalid_argument([&default_not_allowed_schema] {
        methods::validate_method_option_schema(default_not_allowed_schema);
    }));

    return 0;
}
