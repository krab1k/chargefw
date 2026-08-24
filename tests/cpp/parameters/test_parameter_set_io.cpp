#include <chargefw/methods/method_registry.h>
#include <chargefw/parameters/io/parameter_set_io.h>

#include <filesystem>
#include <sstream>
#include <stdexcept>
#include <string>
#include <string_view>
#include <unordered_set>

#include <snitch/snitch.hpp>

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

TEST_CASE("parameter set loads from a PEOE JSON file", "[parameters][io]") {
    const auto parameter_set = parameters::load_parameter_set_json_file(peoe_file());

    CHECK(parameter_set.id() == std::string_view{"PEOE_original"});
    CHECK(parameter_set.method_id() == std::string_view{"peoe"});
    CHECK(parameter_set.name() == std::string_view{"Gasteiger 1980"});
    CHECK(parameter_set.publication() == std::string_view{"10.1016/0040-4020(80)80168-2"});
    CHECK(parameter_set.notes() == std::string_view{"Derived directly from IPs and EAs"});

    CHECK(parameter_set.common().size() == 1);
    CHECK(parameter_set.common().contains("dampH"));
    CHECK(parameter_set.common().parameter("dampH") == 20.02);

    CHECK(parameter_set.atom().size() == 14);

    const auto& hydrogen = parameter_set.atom()[0];

    CHECK(hydrogen.key.atomic_number == 1);
    CHECK(hydrogen.key.classification ==
          parameters::AtomParameterClassificationKind::HIGHEST_BOND_ORDER);
    CHECK(hydrogen.key.type == "1");

    CHECK(parameter_set.atom().contains(0, "A"));
    CHECK(parameter_set.atom().contains(0, "B"));
    CHECK(parameter_set.atom().contains(0, "C"));

    CHECK(parameter_set.atom().parameter(0, "A") == 7.17);
    CHECK(parameter_set.atom().parameter(0, "B") == 6.24);
    CHECK(parameter_set.atom().parameter(0, "C") == -0.56);
}

TEST_CASE("parameter set directory loads all JSON files", "[parameters][io]") {
    const auto parameter_sets =
        parameters::load_parameter_sets_json_directory(CHARGEFW_TEST_PARAMETER_DIR);

    CHECK_FALSE(parameter_sets.empty());
}

TEST_CASE("bundled parameter sets identify registered methods and satisfy their requirements",
          "[parameters][io]") {
    const auto parameter_sets =
        parameters::load_parameter_sets_json_directory(CHARGEFW_TEST_PARAMETER_DIR);
    const auto& registry = chargefw::methods::method_registry();
    auto ids = std::unordered_set<std::string>{};

    for (const auto& parameter_set : parameter_sets) {
        CAPTURE(parameter_set.id());
        CHECK_FALSE(parameter_set.id().empty());
        CHECK(ids.insert(std::string{parameter_set.id()}).second);

        const auto* method = registry.find(parameter_set.method_id());
        REQUIRE(method != nullptr);
        const auto requirements = method->requirements();

        for (const auto name : requirements.common_parameters) {
            CHECK(parameter_set.common().contains(name));
        }

        for (std::size_t entry_index = 0; entry_index < parameter_set.atom().size();
             ++entry_index) {
            for (const auto name : requirements.atom_parameters) {
                CHECK(parameter_set.atom().contains(entry_index, name));
            }
        }

        for (std::size_t entry_index = 0; entry_index < parameter_set.bond().size();
             ++entry_index) {
            for (const auto name : requirements.bond_parameters) {
                CHECK(parameter_set.bond().contains(entry_index, name));
            }
        }
    }
}

TEST_CASE("parameter set rejects malformed JSON", "[parameters][io]") {
    const auto load_bad_length = [] {
        std::istringstream input{bad_length_json()};
        static_cast<void>(parameters::load_parameter_set_json(input));
    };
    CHECK_THROWS_AS(load_bad_length(), std::invalid_argument);

    const auto load_unknown_element = [] {
        std::istringstream input{unknown_element_json()};
        static_cast<void>(parameters::load_parameter_set_json(input));
    };
    CHECK_THROWS_AS(load_unknown_element(), std::invalid_argument);

    const auto load_unknown_classification = [] {
        std::istringstream input{unknown_classification_json()};
        static_cast<void>(parameters::load_parameter_set_json(input));
    };
    CHECK_THROWS_AS(load_unknown_classification(), std::invalid_argument);
}
