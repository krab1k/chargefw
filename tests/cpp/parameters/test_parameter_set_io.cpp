#include <chargefw/parameters/io/parameter_set_io.h>

#include <cassert>
#include <filesystem>
#include <sstream>
#include <stdexcept>
#include <string_view>

namespace parameters = chargefw::parameters;

namespace {

#ifndef CHARGEFW_TEST_PARAMETER_DIR
#error "CHARGEFW_TEST_PARAMETER_DIR must be defined"
#endif

auto peoe_file() -> std::filesystem::path {
    return std::filesystem::path{CHARGEFW_TEST_PARAMETER_DIR} / "PEOE_original.json";
}

auto bad_length_json() -> std::string {
    return R"json(
{
  "metadata": {
    "id": "bad",
    "name": "Bad",
    "method": "peoe"
  },
  "atom": {
    "names": ["A", "B", "C"],
    "data": [
      {
        "key": ["H", "hbo", "1"],
        "value": [7.17, 6.24]
      }
    ]
  }
}
)json";
}

auto unknown_element_json() -> std::string {
    return R"json(
{
  "metadata": {
    "id": "bad",
    "name": "Bad",
    "method": "peoe"
  },
  "atom": {
    "names": ["A"],
    "data": [
      {
        "key": ["Xx", "hbo", "1"],
        "value": [1.0]
      }
    ]
  }
}
)json";
}

auto unknown_classification_json() -> std::string {
    return R"json(
{
  "metadata": {
    "id": "bad",
    "name": "Bad",
    "method": "peoe"
  },
  "atom": {
    "names": ["A"],
    "data": [
      {
        "key": ["H", "unknown", "1"],
        "value": [1.0]
      }
    ]
  }
}
)json";
}

} // namespace

auto main() -> int {
    const auto parameter_set = parameters::load_parameter_set_json_file(peoe_file());

    assert(parameter_set.id() == std::string_view{"PEOE_original"});
    assert(parameter_set.method_id() == std::string_view{"peoe"});
    assert(parameter_set.name() == std::string_view{"Gasteiger 1980"});
    assert(parameter_set.publication() == std::string_view{"10.1016/0040-4020(80)80168-2"});
    assert(parameter_set.notes() == std::string_view{"Derived directly from IPs and EAs"});

    assert(parameter_set.common().size() == 1);
    assert(parameter_set.common().contains("dampH"));
    assert(parameter_set.common().parameter("dampH") == 20.02);

    assert(parameter_set.atom().size() == 14);

    const auto& hydrogen = parameter_set.atom()[0];

    assert(hydrogen.key.atomic_number == 1);
    assert(hydrogen.key.classification ==
           parameters::AtomParameterClassificationKind::HIGHEST_BOND_ORDER);
    assert(hydrogen.key.type == "1");

    assert(parameter_set.atom().contains(0, "A"));
    assert(parameter_set.atom().contains(0, "B"));
    assert(parameter_set.atom().contains(0, "C"));

    assert(parameter_set.atom().parameter(0, "A") == 7.17);
    assert(parameter_set.atom().parameter(0, "B") == 6.24);
    assert(parameter_set.atom().parameter(0, "C") == -0.56);

    const auto parameter_sets =
        parameters::load_parameter_sets_json_directory(CHARGEFW_TEST_PARAMETER_DIR);

    assert(!parameter_sets.empty());

    bool rejected_bad_length = false;

    try {
        std::istringstream input{bad_length_json()};
        [[maybe_unused]] const auto bad = parameters::load_parameter_set_json(input);
    } catch (const std::invalid_argument&) {
        rejected_bad_length = true;
    }

    assert(rejected_bad_length);

    bool rejected_unknown_element = false;

    try {
        std::istringstream input{unknown_element_json()};
        [[maybe_unused]] const auto bad = parameters::load_parameter_set_json(input);
    } catch (const std::invalid_argument&) {
        rejected_unknown_element = true;
    }

    assert(rejected_unknown_element);

    bool rejected_unknown_classification = false;

    try {
        std::istringstream input{unknown_classification_json()};
        [[maybe_unused]] const auto bad = parameters::load_parameter_set_json(input);
    } catch (const std::invalid_argument&) {
        rejected_unknown_classification = true;
    }

    assert(rejected_unknown_classification);

    return 0;
}