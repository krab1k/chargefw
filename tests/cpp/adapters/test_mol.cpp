#include <chargefw/adapters/native/mol2_input.h>
#include <chargefw/adapters/native/mol_input.h>
#include <chargefw/adapters/native/sdf_input.h>
#include <chargefw/core/bond.h>
#include <exception>
#include <filesystem>
#include <fstream>
#include <optional>
#include <snitch/snitch.hpp>
#include <string>

namespace adapters = chargefw::adapters;
namespace mol = chargefw::adapters::native::mol_input;
namespace mol2 = chargefw::adapters::native::mol2_input;
namespace sdf = chargefw::adapters::native::sdf_input;

namespace {

[[nodiscard]] auto fixture(const std::string_view name) -> std::filesystem::path {
    return std::filesystem::path{CHARGEFW_TEST_SOURCE_DIR} / "tests" / "fixtures" / name;
}

} // namespace

TEST_CASE("adapter contract", "[adapters]") {
    {
        std::ifstream input{fixture("synthetic/mol/v2000/charged_atoms.mol")};
        auto reader = mol::MolReader{input, "charged_v2000.mol"};
        const auto result = reader.next();
        REQUIRE(result.has_value());
        const auto& record = *result;
        CHECK(record.molecule.atom_count() == 2);
        CHECK(record.molecule.bond_count() == 1);
        CHECK(record.molecule.conformer_count() == 1);
        CHECK(record.molecule.atom(0).formal_charge() == 1);
        CHECK(record.molecule.atom(1).formal_charge() == -1);
        CHECK(record.molecule.bond(0).order() == chargefw::core::BondOrder::SINGLE);
        CHECK(record.identity.record_id == "charged-v2000");
        CHECK_FALSE(reader.next().has_value());
    }

    {
        std::ifstream input{fixture("synthetic/mol/v3000/charged_atoms.mol")};
        const auto result = mol::parse_mol(
            input, {.source = "charged_v3000.mol", .record_index = 0, .record_id = {}});
        CHECK(result.molecule.atom_count() == 2);
        CHECK(result.molecule.atom(0).formal_charge() == 1);
        CHECK(result.molecule.atom(1).formal_charge() == -1);
        CHECK(result.molecule.bond(0).first_atom_index() == 0);
        CHECK(result.molecule.bond(0).second_atom_index() == 1);
        CHECK(result.diagnostics.empty());
    }

    {
        std::ifstream input{fixture("synthetic/mol2/aromatic.mol2")};
        auto reader = mol2::Mol2Reader{input, "charged_aromatic.mol2"};
        const auto result = reader.next();
        REQUIRE(result.has_value());
        const auto& mol2_record = *result;
        CHECK(mol2_record.molecule.atom_count() == 3);
        CHECK(mol2_record.molecule.bond_count() == 2);
        CHECK(mol2_record.molecule.atom(0).name() == "N1");
        CHECK(mol2_record.molecule.atom(0).formal_charge() == 0);
        CHECK(mol2_record.molecule.bond(0).order() == chargefw::core::BondOrder::SINGLE);
        CHECK(mol2_record.molecule.bond(1).order() == chargefw::core::BondOrder::DOUBLE);
        CHECK(mol2_record.diagnostics.size() == 1);
        CHECK_FALSE(reader.next().has_value());
    }

    {
        std::ifstream input{fixture("synthetic/mol2/malformed_then_valid.mol2")};
        auto reader = mol2::Mol2Reader{input, "malformed_then_valid.mol2"};
        auto rejected = false;
        try {
            [[maybe_unused]] const auto result = reader.next();
        } catch (const std::exception&) {
            rejected = true;
        }
        CHECK(rejected);
    }

    {
        std::ifstream input{fixture("synthetic/mol2/missing_atom_section_then_valid.mol2")};
        auto reader = mol2::Mol2Reader{input, "missing_atom_section_then_valid.mol2"};
        auto rejected = false;
        try {
            [[maybe_unused]] const auto result = reader.next();
        } catch (const std::exception&) {
            rejected = true;
        }
        CHECK(rejected);
    }

    {
        std::ifstream input{fixture("synthetic/sdf/malformed_then_water.sdf")};
        auto reader = sdf::SdfReader{input, "malformed_then_water.sdf"};
        auto rejected = false;
        try {
            [[maybe_unused]] const auto result = reader.next();
        } catch (const std::exception&) {
            rejected = true;
        }
        CHECK(rejected);
    }

    {
        std::ifstream input{fixture("synthetic/sdf/water.sdf")};
        auto reader = sdf::SdfReader{input, "water.sdf"};
        const auto result = reader.next();
        REQUIRE(result.has_value());
        const auto& water_record = *result;
        CHECK(water_record.molecule.atom_count() == 3);
        CHECK(water_record.molecule.bond_count() == 2);
        CHECK_FALSE(reader.next().has_value());
    }

    {
        std::ifstream input{fixture("synthetic/sdf/many_records.sdf")};
        auto reader = sdf::SdfReader{input, "set500.sdf"};
        std::size_t count = 0;

        while (const auto result = reader.next()) {
            ++count;
        }

        CHECK(count == 500);
    }

    {
        std::ifstream input{fixture("synthetic/sdf/v3000.sdf")};
        auto reader = sdf::SdfReader{input, "set_v3000.sdf"};
        const auto result = reader.next();
        REQUIRE(result.has_value());
        const auto& v3000_record = *result;
        CHECK(v3000_record.molecule.atom_count() == 36);
        CHECK(v3000_record.molecule.bond_count() == 38);
    }

    {
        std::ifstream input{fixture("corpus/sdf/heme/ideal.sdf")};
        auto reader = sdf::SdfReader{input, "HEM_ideal.sdf"};
        const auto result = reader.next();
        REQUIRE(result.has_value());
        const auto& hem_record = *result;
        CHECK(hem_record.molecule.atom_count() == 75);
        CHECK(hem_record.molecule.bond_count() == 80);
        CHECK(hem_record.diagnostics.empty());
        CHECK_FALSE(reader.next().has_value());
    }
}
