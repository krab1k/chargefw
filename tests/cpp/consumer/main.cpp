#include <chargefw/adapters/gemmi/mmcif_input.h>
#include <chargefw/adapters/generated_output.h>
#include <chargefw/adapters/native/json_input.h>
#include <chargefw/parameters/io/parameter_set_io.h>

#include <type_traits>

auto main() -> int {
    static_assert(std::is_class_v<chargefw::adapters::gemmi::mmcif_input::MmcifReader>);
    static_assert(std::is_class_v<chargefw::adapters::native::json_input::JsonReader>);
    static_assert(std::is_enum_v<chargefw::adapters::generated_output::Format>);

    return chargefw::parameters::load_default_parameter_sets().empty() ? 1 : 0;
}
