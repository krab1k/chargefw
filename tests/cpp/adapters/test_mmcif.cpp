#include <cassert>
#include <chargefw/adapters/gemmi/mmcif_input.h>
#include <chargefw/core/bond.h>

#include <sstream>
#include <stdexcept>

namespace mmcif = chargefw::adapters::gemmi::mmcif_input;
namespace gemmi_adapter = chargefw::adapters::gemmi;

auto main() -> int {
    std::istringstream input{R"cif(data_first
loop_
_atom_site.group_PDB
_atom_site.id
_atom_site.type_symbol
_atom_site.label_atom_id
_atom_site.label_alt_id
_atom_site.label_comp_id
_atom_site.label_asym_id
_atom_site.label_seq_id
_atom_site.pdbx_PDB_ins_code
_atom_site.Cartn_x
_atom_site.Cartn_y
_atom_site.Cartn_z
_atom_site.occupancy
_atom_site.B_iso_or_equiv
_atom_site.pdbx_formal_charge
_atom_site.auth_seq_id
_atom_site.auth_comp_id
_atom_site.auth_asym_id
_atom_site.auth_atom_id
_atom_site.pdbx_PDB_model_num
ATOM 1 C CA . ALA A 1 ? 0.0 0.0 0.0 1.0 20.0 0 1 ALA A CA 1
HETATM 2 C C1 . LIG A 2 ? 1.0 0.0 0.0 1.0 20.0 0 2 LIG A C1 1
HETATM 3 O O . HOH A 3 ? 2.0 0.0 0.0 1.0 20.0 0 3 HOH A O 1
ATOM 4 C CA . ALA A 1 ? 0.1 0.0 0.0 1.0 20.0 0 1 ALA A CA 2
HETATM 5 C C1 . LIG A 2 ? 1.1 0.0 0.0 1.0 20.0 0 2 LIG A C1 2
HETATM 6 O O . HOH A 3 ? 2.1 0.0 0.0 1.0 20.0 0 3 HOH A O 2
#
data_second
loop_
_atom_site.group_PDB
_atom_site.id
_atom_site.type_symbol
_atom_site.label_atom_id
_atom_site.label_alt_id
_atom_site.label_comp_id
_atom_site.label_asym_id
_atom_site.label_seq_id
_atom_site.pdbx_PDB_ins_code
_atom_site.Cartn_x
_atom_site.Cartn_y
_atom_site.Cartn_z
_atom_site.occupancy
_atom_site.B_iso_or_equiv
_atom_site.pdbx_formal_charge
_atom_site.auth_seq_id
_atom_site.auth_comp_id
_atom_site.auth_asym_id
_atom_site.auth_atom_id
_atom_site.pdbx_PDB_model_num
ATOM 1 O O . HOH B 1 ? 3.0 0.0 0.0 1.0 20.0 0 1 HOH B O 1
#
)cif"};

    auto reader = mmcif::MmcifReader{input, "structure.cif"};
    const auto first = reader.next();
    assert(first.has_value());
    assert(first->identity.source == "structure.cif");
    assert(first->identity.record_index == 0);
    assert(first->identity.record_id == "first");
    assert(first->molecule.atom_count() == 3);
    assert(first->molecule.bond_count() == 0);
    assert(first->molecule.conformer_count() == 2);
    assert(first->molecule.conformer(1)[0].x == 0.1);

    const auto second = reader.next();
    assert(second.has_value());
    assert(second->identity.record_index == 1);
    assert(second->identity.record_id == "second");
    assert(second->molecule.atom_count() == 1);
    assert(!reader.next().has_value());

    std::istringstream deferred_error_input{R"cif(data_valid
loop_
_atom_site.group_PDB
_atom_site.id
_atom_site.type_symbol
_atom_site.label_atom_id
_atom_site.label_alt_id
_atom_site.label_comp_id
_atom_site.label_asym_id
_atom_site.label_seq_id
_atom_site.pdbx_PDB_ins_code
_atom_site.Cartn_x
_atom_site.Cartn_y
_atom_site.Cartn_z
_atom_site.occupancy
_atom_site.B_iso_or_equiv
_atom_site.pdbx_formal_charge
_atom_site.auth_seq_id
_atom_site.auth_comp_id
_atom_site.auth_asym_id
_atom_site.auth_atom_id
_atom_site.pdbx_PDB_model_num
ATOM 1 C CA . ALA A 1 ? 0.0 0.0 0.0 1.0 20.0 0 1 ALA A CA 1
#
data_invalid
loop_
_atom_site.group_PDB
_atom_site.id
_atom_site.type_symbol
_atom_site.label_atom_id
_atom_site.label_alt_id
_atom_site.label_comp_id
_atom_site.label_asym_id
_atom_site.label_seq_id
_atom_site.pdbx_PDB_ins_code
_atom_site.Cartn_x
_atom_site.Cartn_y
_atom_site.Cartn_z
_atom_site.occupancy
_atom_site.B_iso_or_equiv
_atom_site.pdbx_formal_charge
_atom_site.auth_seq_id
_atom_site.auth_comp_id
_atom_site.auth_asym_id
_atom_site.auth_atom_id
_atom_site.pdbx_PDB_model_num
ATOM 1 Xx CA . ALA A 1 ? 0.0 0.0 0.0 1.0 20.0 0 1 ALA A CA 1
#
)cif"};
    auto deferred_error_reader = mmcif::MmcifReader{deferred_error_input};
    assert(deferred_error_reader.next().has_value());
    bool deferred_error = false;
    try {
        static_cast<void>(deferred_error_reader.next());
    } catch (const std::runtime_error&) {
        deferred_error = true;
    }
    assert(deferred_error);

    std::istringstream filtered_input{R"cif(data_filtered
loop_
_atom_site.group_PDB
_atom_site.id
_atom_site.type_symbol
_atom_site.label_atom_id
_atom_site.label_alt_id
_atom_site.label_comp_id
_atom_site.label_asym_id
_atom_site.label_seq_id
_atom_site.pdbx_PDB_ins_code
_atom_site.Cartn_x
_atom_site.Cartn_y
_atom_site.Cartn_z
_atom_site.occupancy
_atom_site.B_iso_or_equiv
_atom_site.pdbx_formal_charge
_atom_site.auth_seq_id
_atom_site.auth_comp_id
_atom_site.auth_asym_id
_atom_site.auth_atom_id
_atom_site.pdbx_PDB_model_num
ATOM 1 C CA . ALA A 1 ? 0.0 0.0 0.0 1.0 20.0 0 1 ALA A CA 1
HETATM 2 C C1 . LIG A 2 ? 1.0 0.0 0.0 1.0 20.0 0 2 LIG A C1 1
HETATM 3 O O . HOH A 3 ? 2.0 0.0 0.0 1.0 20.0 0 3 HOH A O 1
#
)cif"};
    auto filtered = mmcif::MmcifReader{
        filtered_input, {}, {.selection = gemmi_adapter::RecordSelection::polymers}};
    assert(filtered.next()->molecule.atom_count() == 1);

    const auto strategy_input = R"cif(data_connectivity
loop_
_atom_site.group_PDB
_atom_site.id
_atom_site.type_symbol
_atom_site.label_atom_id
_atom_site.label_alt_id
_atom_site.label_comp_id
_atom_site.label_asym_id
_atom_site.label_seq_id
_atom_site.pdbx_PDB_ins_code
_atom_site.Cartn_x
_atom_site.Cartn_y
_atom_site.Cartn_z
_atom_site.occupancy
_atom_site.B_iso_or_equiv
_atom_site.pdbx_formal_charge
_atom_site.auth_seq_id
_atom_site.auth_comp_id
_atom_site.auth_asym_id
_atom_site.auth_atom_id
_atom_site.pdbx_PDB_model_num
HETATM 1 C C1 . LIG A 1 ? 0.0 0.0 0.0 1.0 20.0 0 1 LIG A C1 1
HETATM 2 O O1 . LIG A 1 ? 1.0 0.0 0.0 1.0 20.0 0 1 LIG A O1 1
HETATM 3 C C2 . LIG B 1 ? 2.0 0.0 0.0 1.0 20.0 0 1 LIG B C2 1
HETATM 4 O O2 . LIG B 1 ? 3.0 0.0 0.0 1.0 20.0 0 1 LIG B O2 1
#
loop_
_chem_comp_bond.comp_id
_chem_comp_bond.atom_id_1
_chem_comp_bond.atom_id_2
_chem_comp_bond.value_order
LIG C1 O1 DOUB
LIG C2 O2 SING
#
loop_
_struct_conn.id
_struct_conn.conn_type_id
_struct_conn.ptnr1_label_asym_id
_struct_conn.ptnr1_label_comp_id
_struct_conn.ptnr1_label_seq_id
_struct_conn.ptnr1_label_atom_id
_struct_conn.ptnr2_label_asym_id
_struct_conn.ptnr2_label_comp_id
_struct_conn.ptnr2_label_seq_id
_struct_conn.ptnr2_label_atom_id
link1 covale A LIG 1 O1 B LIG 1 C2
#
)cif";
    const auto read_strategy = [&](const gemmi_adapter::BondStrategy strategy) {
        std::istringstream strategy_stream{strategy_input};
        auto strategy_reader = mmcif::MmcifReader{strategy_stream, {}, {.bond_strategy = strategy}};
        return strategy_reader.next()->molecule;
    };

    assert(read_strategy(gemmi_adapter::BondStrategy::none).bond_count() == 0);
    const auto explicit_molecule = read_strategy(gemmi_adapter::BondStrategy::explicit_bonds);
    assert(explicit_molecule.bond_count() == 3);
    assert(explicit_molecule.bond(0).order() == chargefw::core::BondOrder::DOUBLE);
    assert(explicit_molecule.bond(1).order() == chargefw::core::BondOrder::SINGLE);
    assert(read_strategy(gemmi_adapter::BondStrategy::hybrid).bond_count() == 3);

    const auto duplicate_input = R"cif(data_duplicate
loop_
_atom_site.group_PDB
_atom_site.id
_atom_site.type_symbol
_atom_site.label_atom_id
_atom_site.label_alt_id
_atom_site.label_comp_id
_atom_site.label_asym_id
_atom_site.label_seq_id
_atom_site.pdbx_PDB_ins_code
_atom_site.Cartn_x
_atom_site.Cartn_y
_atom_site.Cartn_z
_atom_site.occupancy
_atom_site.B_iso_or_equiv
_atom_site.pdbx_formal_charge
_atom_site.auth_seq_id
_atom_site.auth_comp_id
_atom_site.auth_asym_id
_atom_site.auth_atom_id
_atom_site.pdbx_PDB_model_num
ATOM 1 N N . ALA A 1 ? 0.0 0.0 0.0 1.0 20.0 0 1 ALA A N 1
ATOM 2 C CA . ALA A 1 ? 1.0 0.0 0.0 1.0 20.0 0 1 ALA A CA 1
#
loop_
_chem_comp_bond.comp_id
_chem_comp_bond.atom_id_1
_chem_comp_bond.atom_id_2
_chem_comp_bond.value_order
ALA N CA DOUB
#
)cif";
    const auto read_duplicate_bond = [&](const gemmi_adapter::BondStrategy strategy) {
        std::istringstream duplicate_stream{duplicate_input};
        auto reader = mmcif::MmcifReader{duplicate_stream, {}, {.bond_strategy = strategy}};
        const auto molecule = reader.next()->molecule;
        assert(molecule.bond_count() == 1);
        return molecule.bond(0);
    };

    assert(read_duplicate_bond(gemmi_adapter::BondStrategy::templates).order() ==
           chargefw::core::BondOrder::SINGLE);
    assert(read_duplicate_bond(gemmi_adapter::BondStrategy::explicit_bonds).order() ==
           chargefw::core::BondOrder::DOUBLE);
    assert(read_duplicate_bond(gemmi_adapter::BondStrategy::hybrid).order() ==
           chargefw::core::BondOrder::DOUBLE);

    return 0;
}
