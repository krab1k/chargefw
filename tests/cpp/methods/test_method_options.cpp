#include <chargefw/methods/method_options.h>

#include <cassert>
#include <stdexcept>
#include <string>
#include <string_view>

namespace methods = chargefw::methods;

auto main() -> int
{
    const methods::MethodOptionSpec schema[]{
        {
            .id = "enabled",
            .description = "Enable the method",
            .type = methods::MethodOptionType::boolean,
            .default_value = true,
            .choices = {}
        },
        {
            .id = "limit",
            .description = "Iteration limit",
            .type = methods::MethodOptionType::integer,
            .default_value = 25,
            .choices = {}
        },
        {
            .id = "radius",
            .description = "Search radius",
            .type = methods::MethodOptionType::floating_point,
            .default_value = 1.5,
            .choices = {}
        },
        {
            .id = "label",
            .description = "Label",
            .type = methods::MethodOptionType::string,
            .default_value = std::string{"default"},
            .choices = {}
        }
    };

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

    bool rejected_missing_option = false;

    try {
        methods::validate_method_options(schema, methods::MethodOptions{});
    } catch (const std::invalid_argument&) {
        rejected_missing_option = true;
    }

    assert(rejected_missing_option);

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

    return 0;
}
