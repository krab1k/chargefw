#include <chargefw/adapters/native/json_input.h>
#include <chargefw/core/bond.h>
#include <snitch/snitch.hpp>

#include <exception>
#include <sstream>
#include <string>

namespace json = chargefw::adapters::native::json_input;

TEST_CASE("JSON input preserves identity, graph, and conformer mapping", "[adapters][json]") {
    {
        std::istringstream input{R"json(
{
  "schema_version": "1.0",
  "molecules": [
    {
      "id": "water-1",
      "name": "water",
      "atoms": [
        { "atomic_number": 8, "formal_charge": 0 },
        { "atomic_number": 1, "formal_charge": 0 },
        { "atomic_number": 1, "formal_charge": 0 }
      ],
      "bonds": [
        { "atoms": [0, 1], "order": 1 },
        { "atoms": [0, 2], "order": 1 }
      ],
      "conformers": [
        {
          "id": "model-1",
          "coordinates": [
            [0.0, 0.0, 0.0],
            [0.9572, 0.0, 0.0],
            [-0.239987, 0.927297, 0.0]
          ]
        }
      ]
    }
  ]
}
)json"};
        auto reader = json::JsonReader{input, "water.json"};
        const auto result = reader.next();
        REQUIRE(result.has_value());
        const auto& record = *result;
        CHECK(record.identity.source == "water.json");
        CHECK(record.identity.record_index == 0);
        CHECK(record.identity.record_id == "water-1");
        CHECK(record.molecule.name() == "water");
        CHECK(record.molecule.atom_count() == 3);
        CHECK(record.molecule.atom(0).atomic_number() == 8);
        CHECK(record.molecule.atom(0).formal_charge() == 0);
        CHECK(record.molecule.atom(0).name().empty());
        CHECK(record.molecule.bond_count() == 2);
        CHECK(record.molecule.bond(0).order() == chargefw::core::BondOrder::SINGLE);
        CHECK(record.molecule.conformer_count() == 1);
        CHECK(record.molecule.conformer(0).name() == "model-1");
        CHECK(record.molecule.conformer(0)[1].x == 0.9572);
        CHECK_FALSE(reader.next().has_value());
    }

    {
        std::istringstream input{R"json(
{
  "schema_version": "1.0",
  "molecules": [
    {
      "id": "invalid",
      "atoms": [{ "atomic_number": 6, "formal_charge": 0 }],
      "bonds": [{ "atoms": [0, 1], "order": 1 }]
    },
    {
      "id": "valid",
      "atoms": [{ "atomic_number": 1, "formal_charge": 0 }]
    }
  ]
}
)json"};
        auto reader = json::JsonReader{input};
        auto rejected = false;
        try {
            [[maybe_unused]] const auto result = reader.next();
        } catch (const std::exception&) {
            rejected = true;
        }
        CHECK(rejected);
    }

    {
        std::istringstream input{R"json(
{
  "schema_version": "2.0",
  "molecules": []
}
)json"};
        auto rejected = false;
        try {
            [[maybe_unused]] auto reader = json::JsonReader{input};
        } catch (const std::runtime_error&) {
            rejected = true;
        }
        CHECK(rejected);
    }
}
