#include <cassert>
#include <chargefw/adapters/gemmi/pdb_input.h>
#include <chargefw/core/bond.h>

#include <sstream>
#include <string>

namespace gemmi_adapter = chargefw::adapters::gemmi;
namespace pdb = gemmi_adapter::pdb_input;

auto main() -> int {
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
    assert(first.has_value());
    assert(first->identity.source == "water.pdb");
    assert(first->identity.record_index == 0);
    assert(!first->identity.record_id.empty());
    assert(first->molecule.atom_count() == 2);
    assert(first->molecule.bond_count() == 0);
    assert(first->molecule.atom(0).atomic_number() == 8);
    assert(first->molecule.atom(1).atomic_number() == 1);
    assert(first->molecule.atom(1).name() == "H1");
    assert(first->molecule.conformer_count() == 2);
    assert(first->molecule.conformer(0)[1].x == 0.957);
    assert(first->mapping.atom_indices[1] == 1);
    assert(first->mapping.conformer_indices[0] == 0);

    assert(first->molecule.conformer(1).name() == "2");
    assert(first->molecule.conformer(1)[0].x == 0.1);
    assert(first->molecule.conformer(1)[1].x == 1.057);
    assert(!reader.next().has_value());

    {
        std::istringstream selection_input{
            R"pdb(ATOM      1  CA  ALA A   1       0.000   0.000   0.000  1.00 20.00           C  
HETATM    2  C1  LIG A   2       1.000   0.000   0.000  1.00 20.00           C  
HETATM    3  O   HOH A   3       2.000   0.000   0.000  1.00 20.00           O  
END
)pdb"};
        auto all_reader = pdb::PdbReader{selection_input};
        assert(all_reader.next()->molecule.atom_count() == 3);

        std::istringstream ligands_input{
            R"pdb(ATOM      1  CA  ALA A   1       0.000   0.000   0.000  1.00 20.00           C  
HETATM    2  C1  LIG A   2       1.000   0.000   0.000  1.00 20.00           C  
HETATM    3  O   HOH A   3       2.000   0.000   0.000  1.00 20.00           O  
END
)pdb"};
        auto ligands_reader =
            pdb::PdbReader{ligands_input, {}, gemmi_adapter::RecordSelection::polymers_and_ligands};
        assert(ligands_reader.next()->molecule.atom_count() == 2);

        std::istringstream polymers_input{
            R"pdb(ATOM      1  CA  ALA A   1       0.000   0.000   0.000  1.00 20.00           C  
HETATM    2  C1  LIG A   2       1.000   0.000   0.000  1.00 20.00           C  
HETATM    3  O   HOH A   3       2.000   0.000   0.000  1.00 20.00           O  
END
)pdb"};
        auto polymers_reader =
            pdb::PdbReader{polymers_input, {}, gemmi_adapter::RecordSelection::polymers};
        assert(polymers_reader.next()->molecule.atom_count() == 1);
    }

    {
        std::istringstream template_input{
            R"pdb(ATOM      1  N   ALA A   1       0.000   0.000   0.000  1.00 20.00           N  
ATOM      2  CA  ALA A   1       1.450   0.000   0.000  1.00 20.00           C  
ATOM      3  C   ALA A   1       2.450   1.000   0.000  1.00 20.00           C  
ATOM      4  O   ALA A   1       3.450   1.000   0.000  1.00 20.00           O  
ATOM      5  CB  ALA A   1       1.450  -1.000   0.000  1.00 20.00           C  
ATOM      6  N   GLY A   2       2.200   2.200   0.000  1.00 20.00           N  
ATOM      7  CA  GLY A   2       3.200   3.200   0.000  1.00 20.00           C  
ATOM      8  C   GLY A   2       4.200   3.200   0.000  1.00 20.00           C  
ATOM      9  O   GLY A   2       5.200   3.200   0.000  1.00 20.00           O  
END
)pdb"};
        auto reader = pdb::PdbReader{template_input,
                                     {},
                                     gemmi_adapter::RecordSelection::all,
                                     gemmi_adapter::BondStrategy::templates};
        const auto record = reader.next();
        assert(record->molecule.bond_count() == 8);
        assert(record->molecule.bond(3).order() == chargefw::core::BondOrder::DOUBLE);
    }

    return 0;
}
