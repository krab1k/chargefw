#include <cassert>
#include <chargefw/adapters/gemmi/mmcif_input.h>

#include <sstream>

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
    auto filtered =
        mmcif::MmcifReader{filtered_input, {}, gemmi_adapter::RecordSelection::polymers};
    assert(filtered.next()->molecule.atom_count() == 1);

    return 0;
}
