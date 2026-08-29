#include <chargefw/adapters/conformer_selection.h>
#include <chargefw/adapters/gemmi/input_options.h>
#include <chargefw/adapters/gemmi/pdb_input.h>
#include <chargefw/core/bond.h>
#include <snitch/snitch.hpp>

#include <exception>
#include <sstream>
#include <string>

namespace gemmi_adapter = chargefw::adapters::gemmi;
namespace pdb = gemmi_adapter::pdb_input;

TEST_CASE("structural input options share stable string conversion", "[adapters][options]") {
    CHECK(gemmi_adapter::record_selection_from_string("polymers-and-ligands") ==
          gemmi_adapter::RecordSelection::polymers_and_ligands);
    CHECK(gemmi_adapter::to_string(gemmi_adapter::RecordSelection::polymers_and_ligands) ==
          "polymers-and-ligands");
    CHECK(gemmi_adapter::bond_strategy_from_string("explicit") ==
          gemmi_adapter::BondStrategy::explicit_bonds);
    CHECK(gemmi_adapter::to_string(gemmi_adapter::BondStrategy::hybrid) == "hybrid");
    CHECK(chargefw::adapters::conformer_selection_from_string("first") ==
          chargefw::adapters::ConformerSelection::first);
    CHECK(chargefw::adapters::to_string(chargefw::adapters::ConformerSelection::all) == "all");

    CHECK_THROWS_AS(gemmi_adapter::record_selection_from_string("unknown"), std::invalid_argument);
}

TEST_CASE("PDB input preserves selected model identity and connectivity", "[adapters][pdb]") {
    std::istringstream input{R"pdb(HEADER    TEST PDB
TITLE     TWO MODELS
MODEL        1
ATOM      1  O   HOH A   1       0.000   0.000   0.000  1.00 20.00           O  
ATOM      2  H1 AHOH A   1       0.957   0.000   0.000  1.00 20.00           H  
ATOM      3  H1 BHOH A   1       9.000   0.000   0.000  1.00 20.00           H  
ENDMDL
MODEL        2
ATOM      1  O   HOH A   1       0.100   0.000   0.000  1.00 20.00           O  
ATOM      2  H1 AHOH A   1       1.057   0.000   0.000  1.00 20.00           H  
ATOM      3  H1 BHOH A   1       9.100   0.000   0.000  1.00 20.00           H  
ENDMDL
END
)pdb"};

    auto reader = pdb::PdbReader{input, "water.pdb"};
    const auto first = reader.next();
    REQUIRE(first.has_value());
    CHECK(first->identity.source == "water.pdb");
    CHECK(first->identity.record_index == 0);
    CHECK_FALSE(first->identity.record_id.empty());
    CHECK(first->molecule.atom_count() == 2);
    CHECK(first->molecule.bond_count() == 0);
    CHECK(first->molecule.atom(0).atomic_number() == 8);
    CHECK(first->molecule.atom(1).atomic_number() == 1);
    CHECK(first->molecule.atom(1).name() == "H1");
    CHECK(first->molecule.conformer_count() == 2);
    CHECK(first->molecule.conformer(0)[1].x == 0.957);

    CHECK(first->molecule.conformer(1).name() == "2");
    CHECK(first->molecule.conformer(1)[0].x == 0.1);
    CHECK(first->molecule.conformer(1)[1].x == 1.057);
    CHECK_FALSE(reader.next().has_value());

    {
        std::istringstream selection_input{
            R"pdb(ATOM      1  CA  ALA A   1       0.000   0.000   0.000  1.00 20.00           C  
HETATM    2  C1  LIG A   2       1.000   0.000   0.000  1.00 20.00           C  
HETATM    3  O   HOH A   3       2.000   0.000   0.000  1.00 20.00           O  
END
)pdb"};
        auto all_reader = pdb::PdbReader{selection_input};
        const auto all_record = all_reader.next();
        REQUIRE(all_record.has_value());
        CHECK(all_record->molecule.atom_count() == 3);

        std::istringstream ligands_input{
            R"pdb(ATOM      1  CA  ALA A   1       0.000   0.000   0.000  1.00 20.00           C  
HETATM    2  C1  LIG A   2       1.000   0.000   0.000  1.00 20.00           C  
HETATM    3  O   HOH A   3       2.000   0.000   0.000  1.00 20.00           O  
END
)pdb"};
        auto ligands_reader = pdb::PdbReader{
            ligands_input, {}, {.selection = gemmi_adapter::RecordSelection::polymers_and_ligands}};
        const auto ligands_record = ligands_reader.next();
        REQUIRE(ligands_record.has_value());
        CHECK(ligands_record->molecule.atom_count() == 2);

        std::istringstream polymers_input{
            R"pdb(ATOM      1  CA  ALA A   1       0.000   0.000   0.000  1.00 20.00           C  
HETATM    2  C1  LIG A   2       1.000   0.000   0.000  1.00 20.00           C  
HETATM    3  O   HOH A   3       2.000   0.000   0.000  1.00 20.00           O  
END
)pdb"};
        auto polymers_reader = pdb::PdbReader{
            polymers_input, {}, {.selection = gemmi_adapter::RecordSelection::polymers}};
        const auto polymers_record = polymers_reader.next();
        REQUIRE(polymers_record.has_value());
        CHECK(polymers_record->molecule.atom_count() == 1);
    }

    {
        const auto strategy_input =
            R"pdb(SSBOND   1 CYS A   3    CYS A   4                          
LINK         C   ALA A   1                 C1  LIG A   5
ATOM      1  N   ALA A   1       0.000   0.000   0.000  1.00 20.00           N  
ATOM      2  CA  ALA A   1       1.450   0.000   0.000  1.00 20.00           C  
ATOM      3  C   ALA A   1       2.450   1.000   0.000  1.00 20.00           C  
ATOM      4  O   ALA A   1       3.450   1.000   0.000  1.00 20.00           O  
ATOM      5  CB  ALA A   1       1.450  -1.000   0.000  1.00 20.00           C  
ATOM      6  N   GLY A   2       2.200   2.200   0.000  1.00 20.00           N  
ATOM      7  CA  GLY A   2       3.200   3.200   0.000  1.00 20.00           C  
ATOM      8  C   GLY A   2       4.200   3.200   0.000  1.00 20.00           C  
ATOM      9  O   GLY A   2       5.200   3.200   0.000  1.00 20.00           O  
ATOM     10  SG  CYS A   3       6.200   3.200   0.000  1.00 20.00           S  
ATOM     11  SG  CYS A   4       7.200   3.200   0.000  1.00 20.00           S  
HETATM   12  C1  LIG A   5       8.200   3.200   0.000  1.00 20.00           C  
CONECT   11   12
END
)pdb";
        const auto read_bond_count = [&](const gemmi_adapter::BondStrategy strategy) {
            std::istringstream strategy_stream{strategy_input};
            auto reader = pdb::PdbReader{strategy_stream, {}, {.bond_strategy = strategy}};
            const auto record = reader.next();
            REQUIRE(record.has_value());
            return record->molecule.bond_count();
        };

        CHECK(read_bond_count(gemmi_adapter::BondStrategy::none) == 0);
        CHECK(read_bond_count(gemmi_adapter::BondStrategy::templates) == 8);
        CHECK(read_bond_count(gemmi_adapter::BondStrategy::explicit_bonds) == 2);
        CHECK(read_bond_count(gemmi_adapter::BondStrategy::hybrid) == 10);
    }

    {
        const auto duplicate_input =
            R"pdb(ATOM      1  N   ALA A   1       0.000   0.000   0.000  1.00 20.00           N
ATOM      2  CA  ALA A   1       1.450   0.000   0.000  1.00 20.00           C
CONECT    1    2
END
)pdb";
        const auto read_bond_count = [&](const gemmi_adapter::BondStrategy strategy) {
            std::istringstream duplicate_stream{duplicate_input};
            auto reader = pdb::PdbReader{duplicate_stream, {}, {.bond_strategy = strategy}};
            const auto record = reader.next();
            REQUIRE(record.has_value());
            return record->molecule.bond_count();
        };

        CHECK(read_bond_count(gemmi_adapter::BondStrategy::templates) == 1);
        CHECK(read_bond_count(gemmi_adapter::BondStrategy::explicit_bonds) == 1);
        CHECK(read_bond_count(gemmi_adapter::BondStrategy::hybrid) == 1);
    }
}

TEST_CASE("PDB input rejects empty and incompatible selected models", "[adapters][pdb]") {
    {
        std::istringstream input{"END\n"};
        bool rejected = false;
        try {
            [[maybe_unused]] auto reader = pdb::PdbReader{input};
        } catch (const std::exception&) {
            rejected = true;
        }
        CHECK(rejected);
    }

    {
        std::istringstream input{R"pdb(MODEL        1
ATOM      1  C   UNL A   1       0.000   0.000   0.000  1.00 20.00           C
ENDMDL
MODEL        2
ATOM      1  O   UNL A   1       0.000   0.000   0.000  1.00 20.00           O
ENDMDL
END
)pdb"};
        bool rejected = false;
        try {
            [[maybe_unused]] auto reader = pdb::PdbReader{input};
        } catch (const std::exception&) {
            rejected = true;
        }
        CHECK(rejected);
    }
}

TEST_CASE("PDB template links respect chains and insertion codes", "[adapters][pdb]") {
    const auto read_bond_count = [](const std::string_view text) {
        std::istringstream input{std::string{text}};
        auto reader =
            pdb::PdbReader{input, {}, {.bond_strategy = gemmi_adapter::BondStrategy::templates}};
        const auto record = reader.next();
        REQUIRE(record.has_value());
        return record->molecule.bond_count();
    };

    CHECK(read_bond_count(
              R"pdb(ATOM      1  C   ALA A   1       0.000   0.000   0.000  1.00 20.00           C
ATOM      2  N   GLY B   2       1.000   0.000   0.000  1.00 20.00           N
END
)pdb") == 0);
    CHECK(read_bond_count(
              R"pdb(ATOM      1  C   ALA A   1       0.000   0.000   0.000  1.00 20.00           C
ATOM      2  N   GLY A   1A      1.000   0.000   0.000  1.00 20.00           N
END
)pdb") == 1);
}
