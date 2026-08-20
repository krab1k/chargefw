#include <cassert>
#include <chargefw/adapters/gemmi/mmcif_input.h>
#include <chargefw/adapters/gemmi/mmcif_output.h>
#include <chargefw/adapters/gemmi/pdb_input.h>
#include <chargefw/charges/atomic_charges.h>
#include <chargefw/charges/charge_collection.h>
#include <chargefw/core/atom.h>
#include <chargefw/core/bond.h>
#include <chargefw/core/conformer.h>
#include <chargefw/core/molecule.h>

#include <gemmi/cif.hpp>
#include <gemmi/read_cif.hpp>

#include <sstream>
#include <vector>

namespace adapters = chargefw::adapters;
namespace charges = chargefw::charges;
namespace core = chargefw::core;
namespace mmcif_input = adapters::gemmi::mmcif_input;
namespace mmcif_output = adapters::gemmi::mmcif_output;

namespace {

auto generated_record(std::string id, const double shift = 0.0)
    -> adapters::ImportedMoleculeRecord {
    auto molecule = core::Molecule{
        {core::Atom{6, 0, "C1"}, core::Atom{8, -1, "O1"}},
        {core::Bond{0, 1, core::BondOrder::DOUBLE}},
        {core::Conformer{{{.x = shift, .y = 0.0, .z = 0.0}, {.x = shift + 1.2, .y = 0.0, .z = 0.0}},
                         "1"}},
        id};
    return {.molecule = std::move(molecule),
            .identity = {.source = "input.sdf", .record_index = 0, .record_id = std::move(id)},
            .diagnostics = {}};
}

auto charge_set(const std::size_t molecule_count) -> charges::ChargeSet {
    std::vector<charges::ChargeAssignment> assignments;
    assignments.reserve(molecule_count);
    for (std::size_t index = 0; index < molecule_count; ++index) {
        assignments.push_back({.target = {.molecule_index = index, .conformer_index = 0},
                               .charges = charges::AtomicCharges{{0.25, -0.25}}});
    }
    return charges::ChargeSet{"qeq", std::move(assignments), "qeq-default"};
}

auto two_conformer_charges() -> charges::ChargeSet {
    return charges::ChargeSet{"qeq",
                              {{.target = {.molecule_index = 0, .conformer_index = 0},
                                .charges = charges::AtomicCharges{{0.25, -0.25}}},
                               {.target = {.molecule_index = 0, .conformer_index = 1},
                                .charges = charges::AtomicCharges{{0.20, -0.20}}}},
                              "qeq-default"};
}

} // namespace

auto main() -> int {
    {
        auto records = std::vector<adapters::ImportedMoleculeRecord>{generated_record("first"),
                                                                     generated_record("second", 2)};
        std::ostringstream output;
        mmcif_output::MmcifWriter{output}.write_generated(records, charge_set(records.size()));
        const auto document = ::gemmi::cif::read_string(output.str());
        assert(document.blocks.size() == 2);
        assert(document.blocks[0].has_mmcif_category("_chem_comp."));
        assert(document.blocks[1].has_mmcif_category("_chem_comp_bond."));
        assert(document.blocks[0].has_mmcif_category("_sb_ncbr_partial_atomic_charges_meta."));

        std::istringstream round_trip{output.str()};
        auto reader = mmcif_input::MmcifReader{
            round_trip, {}, {.bond_strategy = adapters::gemmi::BondStrategy::explicit_bonds}};
        const auto first = reader.next();
        const auto second = reader.next();
        assert(first->molecule.atom_count() == 2);
        assert(first->molecule.bond_count() == 1);
        assert(first->molecule.bond(0).order() == core::BondOrder::DOUBLE);
        assert(first->molecule.atom(1).formal_charge() == -1);
        assert(second->molecule.conformer(0)[0].x == 2.0);
    }

    {
        std::istringstream input{R"cif(data_source
_entry.id source
_custom_extension.note 'keep me'
_audit_conform.dict_name mmcif_pdbx.dic
_audit_conform.dict_version 5.0
_audit_conform.dict_location https://mmcif.wwpdb.org
loop_
_atom_site.group_PDB
_atom_site.id
_atom_site.type_symbol
_atom_site.label_atom_id
_atom_site.label_alt_id
_atom_site.label_comp_id
_atom_site.label_asym_id
_atom_site.label_seq_id
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
HETATM 1 C C1 . UNL A 1 0 0 0 1 0 0 1 UNL A C1 1
HETATM 2 O O1 . UNL A 1 1 0 0 1 0 -1 1 UNL A O1 1
#
)cif"};
        auto reader = mmcif_input::MmcifReader{input};
        auto records = std::vector<adapters::ImportedMoleculeRecord>{std::move(*reader.next())};
        const auto source =
            mmcif_output::MmcifSource{.document = reader.source_document(),
                                      .block_indices = reader.source_block_indices(),
                                      .selection = reader.options().selection};
        std::ostringstream output;
        mmcif_output::MmcifWriter{output}.write_mmcif(records, charge_set(1), source);
        const auto document = ::gemmi::cif::read_string(output.str());
        assert(::gemmi::cif::as_string(*document.blocks[0].find_value("_custom_extension.note")) ==
               "keep me");
        assert(document.blocks[0].has_mmcif_category("_sb_ncbr_partial_atomic_charges."));

        std::istringstream appended_input{output.str()};
        auto appended_reader = mmcif_input::MmcifReader{appended_input};
        auto appended_records =
            std::vector<adapters::ImportedMoleculeRecord>{std::move(*appended_reader.next())};
        const auto appended_source =
            mmcif_output::MmcifSource{.document = appended_reader.source_document(),
                                      .block_indices = appended_reader.source_block_indices(),
                                      .selection = appended_reader.options().selection};
        std::ostringstream appended_output;
        mmcif_output::MmcifWriter{appended_output}.write_mmcif(appended_records, charge_set(1),
                                                               appended_source, "ChargeFW", "0.0.1",
                                                               mmcif_output::WriteMode::append);
        auto appended = ::gemmi::cif::read_string(appended_output.str());
        assert(appended.blocks[0]
                   .find_mmcif_category("_sb_ncbr_partial_atomic_charges_meta.")
                   .length() == 2);

        std::istringstream replaced_input{appended_output.str()};
        auto replaced_reader = mmcif_input::MmcifReader{replaced_input};
        auto replaced_records =
            std::vector<adapters::ImportedMoleculeRecord>{std::move(*replaced_reader.next())};
        const auto replaced_source =
            mmcif_output::MmcifSource{.document = replaced_reader.source_document(),
                                      .block_indices = replaced_reader.source_block_indices(),
                                      .selection = replaced_reader.options().selection};
        std::ostringstream replaced_output;
        mmcif_output::MmcifWriter{replaced_output}.write_mmcif(replaced_records, charge_set(1),
                                                               replaced_source);
        auto replaced = ::gemmi::cif::read_string(replaced_output.str());
        assert(replaced.blocks[0]
                   .find_mmcif_category("_sb_ncbr_partial_atomic_charges_meta.")
                   .length() == 1);
    }

    {
        std::istringstream input{R"cif(data_filtered
loop_
_atom_site.group_PDB
_atom_site.id
_atom_site.type_symbol
_atom_site.label_atom_id
_atom_site.label_alt_id
_atom_site.label_comp_id
_atom_site.label_asym_id
_atom_site.label_seq_id
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
ATOM 1 C CA . ALA A 1 0 0 0 1 0 0 1 ALA A CA 1
HETATM 2 C C1 . LIG A 2 1 0 0 1 0 0 2 LIG A C1 1
HETATM 3 O O . HOH A 3 2 0 0 1 0 0 3 HOH A O 1
#
)cif"};
        auto reader = mmcif_input::MmcifReader{
            input, {}, {.selection = adapters::gemmi::RecordSelection::polymers}};
        auto records = std::vector<adapters::ImportedMoleculeRecord>{std::move(*reader.next())};
        const auto source =
            mmcif_output::MmcifSource{.document = reader.source_document(),
                                      .block_indices = reader.source_block_indices(),
                                      .selection = reader.options().selection};
        const auto charges = charges::ChargeSet{
            "formal",
            {{.target = {.molecule_index = 0}, .charges = charges::AtomicCharges{{0.0}}}}};
        std::ostringstream output;
        mmcif_output::MmcifWriter{output}.write_mmcif(records, charges, source);
        auto document = ::gemmi::cif::read_string(output.str());
        auto charge_rows =
            document.blocks[0].find("_sb_ncbr_partial_atomic_charges.", {"atom_id", "charge"});
        assert(charge_rows.length() == 1);
        assert(::gemmi::cif::as_string(charge_rows[0][0]) == "1");
    }

    {
        std::istringstream input{R"pdb(MODEL        1
HETATM    1  C1  UNL A   1       0.000   0.000   0.000  1.00  0.00           C
HETATM    2  O1  UNL A   1       1.200   0.000   0.000  1.00  0.00           O1-
ENDMDL
MODEL        2
HETATM    1  C1  UNL A   1       0.100   0.000   0.000  1.00  0.00           C
HETATM    2  O1  UNL A   1       1.300   0.000   0.000  1.00  0.00           O1-
ENDMDL
END
)pdb"};
        auto reader = adapters::gemmi::pdb_input::PdbReader{input, "models.pdb"};
        auto records = std::vector<adapters::ImportedMoleculeRecord>{std::move(*reader.next())};
        const auto source = mmcif_output::PdbSource{.structure = reader.source_structure(),
                                                    .selection = reader.options().selection};
        std::ostringstream output;
        mmcif_output::MmcifWriter{output}.write_pdb(records.front(), two_conformer_charges(),
                                                    source);
        auto document = ::gemmi::cif::read_string(output.str());
        auto metadata = document.blocks[0].find("_sb_ncbr_partial_atomic_charges_meta.",
                                                {"id", "type", "method"});
        assert(metadata.length() == 2);
        assert(::gemmi::cif::as_string(metadata[0][1]) == "empirical");
        assert(::gemmi::cif::as_string(metadata[1][2]) == "qeq/qeq-default");

        auto charges = document.blocks[0].find("_sb_ncbr_partial_atomic_charges.",
                                               {"type_id", "atom_id", "charge"});
        assert(charges.length() == 4);
        assert(::gemmi::cif::as_string(charges[0][0]) != ::gemmi::cif::as_string(charges[2][0]));
        assert(::gemmi::cif::as_string(charges[0][1]) != ::gemmi::cif::as_string(charges[2][1]));

        const auto no_parameter_charges =
            charges::ChargeSet{"formal",
                               {{.target = {.molecule_index = 0, .conformer_index = 0},
                                 .charges = charges::AtomicCharges{{0.0, 0.0}}}}};
        std::ostringstream no_parameter_output;
        mmcif_output::MmcifWriter{no_parameter_output}.write_pdb(records.front(),
                                                                 no_parameter_charges, source);
        auto no_parameter_document = ::gemmi::cif::read_string(no_parameter_output.str());
        auto no_parameter_metadata = no_parameter_document.blocks[0].find(
            "_sb_ncbr_partial_atomic_charges_meta.", {"method"});
        assert(::gemmi::cif::as_string(no_parameter_metadata[0][0]) == "formal");

        std::istringstream round_trip{output.str()};
        auto round_trip_reader = mmcif_input::MmcifReader{round_trip};
        assert(round_trip_reader.next()->molecule.conformer_count() == 2);
    }

    return 0;
}
