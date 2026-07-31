#include <cassert>
#include <chargefw/adapters/native/mol2_input.h>
#include <chargefw/adapters/native/mol_input.h>
#include <chargefw/adapters/native/sdf_input.h>
#include <chargefw/core/bond.h>
#include <exception>
#include <filesystem>
#include <fstream>
#include <optional>
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

auto main() -> int {
    {
        std::ifstream input{fixture("synthetic/mol/v2000/charged_atoms.mol")};
        auto reader = mol::MolReader{input, "charged_v2000.mol"};
        const auto result = reader.next();
        if (!result.has_value()) {
            return 1;
        }
        const auto& record = *result;
        assert(record.molecule.atom_count() == 2);
        assert(record.molecule.bond_count() == 1);
        assert(record.molecule.conformer_count() == 1);
        assert(record.molecule.atom(0).formal_charge() == 1);
        assert(record.molecule.atom(1).formal_charge() == -1);
        assert(record.molecule.bond(0).order() == chargefw::core::BondOrder::SINGLE);
        assert(record.identity.record_id == "charged-v2000");
        assert(!reader.next().has_value());
    }

    {
        std::ifstream input{fixture("synthetic/mol/v3000/charged_atoms.mol")};
        const auto result = mol::parse_mol(
            input, {.source = "charged_v3000.mol", .record_index = 0, .record_id = {}});
        assert(result.molecule.atom_count() == 2);
        assert(result.molecule.atom(0).formal_charge() == 1);
        assert(result.molecule.atom(1).formal_charge() == -1);
        assert(result.molecule.bond(0).first_atom_index() == 0);
        assert(result.molecule.bond(0).second_atom_index() == 1);
        assert(result.diagnostics.empty());
    }

    {
        std::ifstream input{fixture("synthetic/mol2/aromatic.mol2")};
        auto reader = mol2::Mol2Reader{input, "charged_aromatic.mol2"};
        const auto result = reader.next();
        if (!result.has_value()) {
            return 1;
        }
        const auto& mol2_record = *result;
        assert(mol2_record.molecule.atom_count() == 3);
        assert(mol2_record.molecule.bond_count() == 2);
        assert(mol2_record.molecule.atom(0).name() == "N1");
        assert(mol2_record.molecule.atom(0).formal_charge() == 0);
        assert(mol2_record.molecule.bond(0).order() == chargefw::core::BondOrder::AROMATIC);
        assert(mol2_record.molecule.bond(1).order() == chargefw::core::BondOrder::DOUBLE);
        assert(mol2_record.diagnostics.size() == 1);
        assert(!reader.next().has_value());
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
        assert(rejected);
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
        assert(rejected);
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
        assert(rejected);
    }

    {
        std::ifstream input{fixture("synthetic/sdf/water.sdf")};
        auto reader = sdf::SdfReader{input, "water.sdf"};
        const auto result = reader.next();
        if (!result.has_value()) {
            return 1;
        }
        const auto& water_record = *result;
        assert(water_record.molecule.atom_count() == 3);
        assert(water_record.molecule.bond_count() == 2);
        assert(!reader.next().has_value());
    }

    {
        std::ifstream input{fixture("synthetic/sdf/many_records.sdf")};
        auto reader = sdf::SdfReader{input, "set500.sdf"};
        std::size_t count = 0;

        while (const auto result = reader.next()) {
            ++count;
        }

        assert(count == 500);
    }

    {
        std::ifstream input{fixture("synthetic/sdf/v3000.sdf")};
        auto reader = sdf::SdfReader{input, "set_v3000.sdf"};
        const auto result = reader.next();
        if (!result.has_value()) {
            return 1;
        }
        const auto& v3000_record = *result;
        assert(v3000_record.molecule.atom_count() == 36);
        assert(v3000_record.molecule.bond_count() == 38);
    }

    {
        std::ifstream input{fixture("corpus/sdf/heme/ideal.sdf")};
        auto reader = sdf::SdfReader{input, "HEM_ideal.sdf"};
        const auto result = reader.next();
        if (!result.has_value()) {
            return 1;
        }
        const auto& hem_record = *result;
        assert(hem_record.molecule.atom_count() == 75);
        assert(hem_record.molecule.bond_count() == 80);
        assert(hem_record.diagnostics.empty());
        assert(!reader.next().has_value());
    }

    return 0;
}
