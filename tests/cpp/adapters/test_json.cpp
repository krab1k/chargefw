#include <cassert>
#include <chargefw/adapters/native/json_input.h>
#include <chargefw/core/bond.h>

#include <sstream>
#include <string>

namespace json = chargefw::adapters::native::json_input;

auto main() -> int {
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
        assert(result.has_value());
        assert(result->has_value());
        const auto& record = result->value();
        assert(record.identity.source == "water.json");
        assert(record.identity.record_index == 0);
        assert(record.identity.record_id == "water-1");
        assert(record.molecule.name() == "water");
        assert(record.molecule.atom_count() == 3);
        assert(record.molecule.atom(0).atomic_number() == 8);
        assert(record.molecule.atom(0).formal_charge() == 0);
        assert(record.molecule.atom(0).name().empty());
        assert(record.molecule.bond_count() == 2);
        assert(record.molecule.bond(0).order() == chargefw::core::BondOrder::SINGLE);
        assert(record.molecule.conformer_count() == 1);
        assert(record.molecule.conformer(0).name() == "model-1");
        assert(record.molecule.conformer(0)[1].x == 0.9572);
        assert(record.mapping.atom_indices[2] == 2);
        assert(record.mapping.conformer_indices[0] == 0);
        assert(!reader.next().has_value());
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
        const auto first = reader.next();
        const auto second = reader.next();
        assert(first.has_value());
        assert(!first->has_value());
        assert(first->error().identity.record_index == 0);
        assert(second.has_value());
        assert(second->has_value());
        assert(second->value().identity.record_id == "valid");
        assert(!reader.next().has_value());
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
        assert(rejected);
    }

    return 0;
}
