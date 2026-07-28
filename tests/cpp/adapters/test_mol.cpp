#include <cassert>
#include <chargefw/adapters/native/mol.h>
#include <chargefw/adapters/native/mol2.h>
#include <chargefw/adapters/native/sdf.h>
#include <chargefw/core/bond.h>
#include <filesystem>
#include <fstream>
#include <optional>
#include <string>

namespace adapters = chargefw::adapters;
namespace mol = chargefw::adapters::native::mol;
namespace mol2 = chargefw::adapters::native::mol2;
namespace sdf = chargefw::adapters::native::sdf;

namespace {

[[nodiscard]] auto fixture(const std::string_view name) -> std::filesystem::path {
    return std::filesystem::path{CHARGEFW_TEST_SOURCE_DIR} / "tests" / "fixtures" / name;
}

} // namespace

auto main() -> int {
    {
        std::ifstream input{fixture("mol/charged_v2000.mol")};
        auto reader = mol::MolReader{input, "charged_v2000.mol"};
        const auto result = reader.next();
        assert(result.has_value());
        assert(result->has_value());
        assert(result->value().molecule.atom_count() == 2);
        assert(result->value().molecule.bond_count() == 1);
        assert(result->value().molecule.conformer_count() == 1);
        assert(result->value().molecule.atom(0).formal_charge() == 1);
        assert(result->value().molecule.atom(1).formal_charge() == -1);
        assert(result->value().molecule.bond(0).order() == chargefw::core::BondOrder::SINGLE);
        assert(result->value().identity.record_id == "charged-v2000");
        assert(!reader.next().has_value());
    }

    {
        std::ifstream input{fixture("mol/charged_v3000.mol")};
        const auto result = mol::parse_mol(
            input, {.source = "charged_v3000.mol", .record_index = 0, .record_id = {}});
        assert(result.has_value());
        assert(result->molecule.atom_count() == 2);
        assert(result->molecule.atom(0).formal_charge() == 1);
        assert(result->molecule.atom(1).formal_charge() == -1);
        assert(result->molecule.bond(0).first_atom_index() == 0);
        assert(result->molecule.bond(0).second_atom_index() == 1);
        assert(result->diagnostics.size() == 1);
    }

    {
        std::ifstream input{fixture("mol2/charged_aromatic.mol2")};
        auto reader = mol2::Mol2Reader{input, "charged_aromatic.mol2"};
        const auto result = reader.next();
        assert(result.has_value());
        assert(result->has_value());
        assert(result->value().molecule.atom_count() == 3);
        assert(result->value().molecule.bond_count() == 2);
        assert(result->value().molecule.atom(0).name() == "N1");
        assert(result->value().molecule.atom(0).formal_charge() == 0);
        assert(result->value().molecule.bond(0).order() == chargefw::core::BondOrder::AROMATIC);
        assert(result->value().molecule.bond(1).order() == chargefw::core::BondOrder::DOUBLE);
        assert(result->value().diagnostics.size() == 1);
        assert(!reader.next().has_value());
    }

    {
        std::ifstream input{fixture("mol2/malformed_then_valid.mol2")};
        auto reader = mol2::Mol2Reader{input, "malformed_then_valid.mol2"};
        const auto first = reader.next();
        const auto second = reader.next();
        assert(first.has_value());
        assert(!first->has_value());
        assert(second.has_value());
        assert(second->has_value());
        assert(second->value().identity.record_index == 1);
        assert(second->value().molecule.atom_count() == 1);
    }

    {
        std::ifstream input{fixture("sdf/malformed_then_water.sdf")};
        auto reader = sdf::SdfReader{input, "malformed_then_water.sdf"};
        const auto first = reader.next();
        const auto second = reader.next();
        assert(first.has_value());
        assert(!first->has_value());
        assert(first->error().identity.record_index == 0);
        assert(second.has_value());
        assert(second->has_value());
        assert(second->value().identity.record_index == 1);
        assert(second->value().molecule.atom_count() == 3);
        assert(!reader.next().has_value());
    }

    {
        std::ifstream input{std::filesystem::path{CHARGEFW_TEST_SOURCE_DIR} / "tests" /
                            "water.sdf"};
        auto reader = sdf::SdfReader{input, "water.sdf"};
        const auto result = reader.next();
        assert(result.has_value());
        assert(result->has_value());
        assert(result->value().molecule.atom_count() == 3);
        assert(result->value().molecule.bond_count() == 2);
        assert(!reader.next().has_value());
    }

    {
        std::ifstream input{std::filesystem::path{CHARGEFW_TEST_SOURCE_DIR} / "tests" /
                            "set500.sdf"};
        auto reader = sdf::SdfReader{input, "set500.sdf"};
        std::size_t count = 0;

        while (const auto result = reader.next()) {
            assert(result->has_value());
            ++count;
        }

        assert(count == 500);
    }

    {
        std::ifstream input{std::filesystem::path{CHARGEFW_TEST_SOURCE_DIR} / "tests" /
                            "set_v3000.sdf"};
        auto reader = sdf::SdfReader{input, "set_v3000.sdf"};
        const auto result = reader.next();
        assert(result.has_value());
        assert(result->has_value());
        assert(result->value().molecule.atom_count() == 36);
        assert(result->value().molecule.bond_count() == 38);
    }

    return 0;
}
