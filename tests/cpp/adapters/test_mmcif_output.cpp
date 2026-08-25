#include <chargefw/adapters/gemmi/mmcif_input.h>
#include <chargefw/adapters/gemmi/mmcif_output.h>
#include <chargefw/adapters/gemmi/pdb_input.h>
#include <chargefw/charges/atomic_charges.h>
#include <chargefw/charges/charge_collection.h>
#include <chargefw/core/atom.h>
#include <chargefw/core/bond.h>
#include <chargefw/core/conformer.h>
#include <chargefw/core/molecule.h>
#include <snitch/snitch.hpp>

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

TEST_CASE("mmCIF generated output preserves charge mapping", "[adapters][mmcif]") {
    auto records = std::vector<adapters::ImportedMoleculeRecord>{generated_record("first"),
                                                                 generated_record("second", 2)};
    std::ostringstream output;
    mmcif_output::MmcifWriter{output}.write_generated(records, charge_set(records.size()));
    auto document = ::gemmi::cif::read_string(output.str());
    REQUIRE(document.blocks.size() == 2);
    CHECK(document.blocks[0].has_mmcif_category("_chem_comp."));
    CHECK(document.blocks[1].has_mmcif_category("_chem_comp_bond."));
    auto atom_types = document.blocks[0].find("_atom_type.", {"symbol"});
    REQUIRE(atom_types.length() == 2);
    CHECK(::gemmi::cif::as_string(atom_types[0][0]) == "C");
    CHECK(::gemmi::cif::as_string(atom_types[1][0]) == "O");
    CHECK(document.blocks[0].has_mmcif_category("_sb_ncbr_partial_atomic_charges_meta."));
    auto dictionary =
        document.blocks[0].find("_audit_conform.", {"dict_name", "dict_version", "dict_location"});
    REQUIRE(dictionary.length() == 1);
    CHECK(::gemmi::cif::as_string(dictionary[0][0]) == "mmcif_charges_v11.dic");
    CHECK(::gemmi::cif::as_string(dictionary[0][1]) == "1.1");
    CHECK(::gemmi::cif::as_string(dictionary[0][2]) ==
          "https://sb-ncbr.github.io/charges-schema/schemas/mmcif_charges_v11.dic");
    auto metadata = document.blocks[0].find(
        "_sb_ncbr_partial_atomic_charges_meta.",
        {"type", "method", "parameter_set", "software_name", "software_version"});
    REQUIRE(metadata.length() == 1);
    CHECK(::gemmi::cif::as_string(metadata[0][0]) == "empirical");
    CHECK(::gemmi::cif::as_string(metadata[0][1]) == "qeq");
    CHECK(::gemmi::cif::as_string(metadata[0][2]) == "qeq-default");
    CHECK(::gemmi::cif::as_string(metadata[0][3]) == "ChargeFW");
    CHECK(::gemmi::cif::as_string(metadata[0][4]) == "unknown");
    const auto* bond_order = document.blocks[0].find_value("_chem_comp_bond.value_order");
    REQUIRE(bond_order != nullptr);
    CHECK(::gemmi::cif::as_string(*bond_order) == "doub");

    std::istringstream round_trip{output.str()};
    auto reader = mmcif_input::MmcifReader{
        round_trip, {}, {.bond_strategy = adapters::gemmi::BondStrategy::explicit_bonds}};
    const auto first = reader.next();
    const auto second = reader.next();
    REQUIRE(first.has_value());
    REQUIRE(second.has_value());
    CHECK(first->molecule.atom_count() == 2);
    CHECK(first->molecule.bond_count() == 1);
    CHECK(first->molecule.bond(0).order() == core::BondOrder::DOUBLE);
    CHECK(first->molecule.atom(1).formal_charge() == -1);
    CHECK(second->molecule.conformer(0)[0].x == 2.0);
}

TEST_CASE("mmCIF output preserves source categories and charge assignments", "[adapters][mmcif]") {
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
    const auto source = mmcif_output::MmcifSource{.document = reader.source_document(),
                                                  .block_indices = reader.source_block_indices(),
                                                  .selection = reader.options().selection};
    std::ostringstream output;
    mmcif_output::MmcifWriter{output}.write_mmcif(records, charge_set(1), source);
    auto document = ::gemmi::cif::read_string(output.str());
    REQUIRE(document.blocks.size() == 1);
    const auto* custom_note = document.blocks[0].find_value("_custom_extension.note");
    REQUIRE(custom_note != nullptr);
    CHECK(::gemmi::cif::as_string(*custom_note) == "keep me");
    CHECK(document.blocks[0].has_mmcif_category("_sb_ncbr_partial_atomic_charges."));

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
    REQUIRE(appended.blocks.size() == 1);
    CHECK(
        appended.blocks[0].find_mmcif_category("_sb_ncbr_partial_atomic_charges_meta.").length() ==
        2);

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
    REQUIRE(replaced.blocks.size() == 1);
    CHECK(
        replaced.blocks[0].find_mmcif_category("_sb_ncbr_partial_atomic_charges_meta.").length() ==
        1);
}

TEST_CASE("mmCIF output maps filtered records to source atom IDs", "[adapters][mmcif]") {
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
    const auto source = mmcif_output::MmcifSource{.document = reader.source_document(),
                                                  .block_indices = reader.source_block_indices(),
                                                  .selection = reader.options().selection};
    const auto charges = charges::ChargeSet{
        "formal", {{.target = {.molecule_index = 0}, .charges = charges::AtomicCharges{{0.0}}}}};
    std::ostringstream output;
    mmcif_output::MmcifWriter{output}.write_mmcif(records, charges, source);
    auto document = ::gemmi::cif::read_string(output.str());
    REQUIRE(document.blocks.size() == 1);
    auto charge_rows =
        document.blocks[0].find("_sb_ncbr_partial_atomic_charges.", {"atom_id", "charge"});
    REQUIRE(charge_rows.length() == 1);
    CHECK(::gemmi::cif::as_string(charge_rows[0][0]) == "1");
}

TEST_CASE("mmCIF output preserves omitted alternate locations", "[adapters][mmcif]") {
    std::istringstream input{R"cif(data_alternate_locations
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
HETATM 11 C C1 B LIG A 1 0 0 0 1 0 0 1 LIG A C1 1
HETATM 12 C C1 A LIG A 1 1 0 0 1 0 0 1 LIG A C1 1
HETATM 13 O O1 . LIG A 1 2 0 0 1 0 -1 1 LIG A O1 1
#
)cif"};
    auto reader = mmcif_input::MmcifReader{input};
    auto records = std::vector<adapters::ImportedMoleculeRecord>{std::move(*reader.next())};
    REQUIRE(records.front().molecule.atom_count() == 2);
    CHECK(records.front().molecule.conformer(0)[0].x == 1.0);
    const auto source = mmcif_output::MmcifSource{.document = reader.source_document(),
                                                  .block_indices = reader.source_block_indices(),
                                                  .selection = reader.options().selection};
    std::ostringstream output;
    mmcif_output::MmcifWriter{output}.write_mmcif(records, charge_set(1), source);
    auto document = ::gemmi::cif::read_string(output.str());
    auto atom_sites =
        document.blocks.front().find("_atom_site.", {"id", "label_atom_id", "label_alt_id"});
    REQUIRE(atom_sites.length() == 3);
    CHECK(::gemmi::cif::as_string(atom_sites[0][0]) == "11");
    CHECK(::gemmi::cif::as_string(atom_sites[0][2]) == "B");
    CHECK(::gemmi::cif::as_string(atom_sites[1][0]) == "12");
    CHECK(::gemmi::cif::as_string(atom_sites[1][2]) == "A");

    auto charge_rows =
        document.blocks.front().find("_sb_ncbr_partial_atomic_charges.", {"atom_id", "charge"});
    REQUIRE(charge_rows.length() == 2);
    CHECK(::gemmi::cif::as_string(charge_rows[0][0]) == "12");
    CHECK(::gemmi::cif::as_string(charge_rows[1][0]) == "13");
}

TEST_CASE("mmCIF output from PDB input preserves selected conformer charge mapping",
          "[adapters][pdb][mmcif]") {
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
    mmcif_output::MmcifWriter{output}.write_pdb(records.front(), two_conformer_charges(), source,
                                                "ChargeFW", "0.0.1");
    auto document = ::gemmi::cif::read_string(output.str());
    REQUIRE(document.blocks.size() == 1);
    auto metadata = document.blocks[0].find(
        "_sb_ncbr_partial_atomic_charges_meta.",
        {"id", "type", "method", "parameter_set", "software_name", "software_version"});
    REQUIRE(metadata.length() == 2);
    CHECK(::gemmi::cif::as_string(metadata[0][1]) == "empirical");
    CHECK(::gemmi::cif::as_string(metadata[1][2]) == "qeq");
    CHECK(::gemmi::cif::as_string(metadata[1][3]) == "qeq-default");
    CHECK(::gemmi::cif::as_string(metadata[1][4]) == "ChargeFW");
    CHECK(::gemmi::cif::as_string(metadata[1][5]) == "0.0.1");

    auto charges = document.blocks[0].find("_sb_ncbr_partial_atomic_charges.",
                                           {"type_id", "atom_id", "charge"});
    REQUIRE(charges.length() == 4);
    CHECK(::gemmi::cif::as_string(charges[0][0]) != ::gemmi::cif::as_string(charges[2][0]));
    CHECK(::gemmi::cif::as_string(charges[0][1]) != ::gemmi::cif::as_string(charges[2][1]));

    const auto no_parameter_charges =
        charges::ChargeSet{"formal",
                           {{.target = {.molecule_index = 0, .conformer_index = 0},
                             .charges = charges::AtomicCharges{{0.0, 0.0}}}}};
    std::ostringstream no_parameter_output;
    mmcif_output::MmcifWriter{no_parameter_output}.write_pdb(records.front(), no_parameter_charges,
                                                             source);
    auto no_parameter_document = ::gemmi::cif::read_string(no_parameter_output.str());
    REQUIRE(no_parameter_document.blocks.size() == 1);
    auto no_parameter_metadata = no_parameter_document.blocks[0].find(
        "_sb_ncbr_partial_atomic_charges_meta.",
        {"method", "parameter_set", "software_name", "software_version"});
    REQUIRE(no_parameter_metadata.length() == 1);
    CHECK(::gemmi::cif::as_string(no_parameter_metadata[0][0]) == "formal");
    CHECK(::gemmi::cif::as_string(no_parameter_metadata[0][1]).empty());
    CHECK(::gemmi::cif::as_string(no_parameter_metadata[0][2]) == "ChargeFW");
    CHECK(::gemmi::cif::as_string(no_parameter_metadata[0][3]) == "unknown");

    std::istringstream round_trip{output.str()};
    auto round_trip_reader = mmcif_input::MmcifReader{round_trip};
    const auto round_trip_record = round_trip_reader.next();
    REQUIRE(round_trip_record.has_value());
    CHECK(round_trip_record->molecule.conformer_count() == 2);
}

TEST_CASE("mmCIF output from PDB input omits unselected alternate locations",
          "[adapters][pdb][mmcif]") {
    std::istringstream input{
        R"pdb(HETATM    1  C1 BLIG A   1       0.000   0.000   0.000  1.00  0.00           C
HETATM    2  C1 ALIG A   1       1.000   0.000   0.000  1.00  0.00           C
HETATM    3  O1  LIG A   1       2.000   0.000   0.000  1.00  0.00           O
END
)pdb"};
    auto reader = adapters::gemmi::pdb_input::PdbReader{input, "alternate-locations.pdb"};
    auto record = std::move(*reader.next());
    REQUIRE(record.molecule.atom_count() == 2);
    CHECK(record.molecule.conformer(0)[0].x == 1.0);
    const auto source = mmcif_output::PdbSource{.structure = reader.source_structure(),
                                                .selection = reader.options().selection};
    std::ostringstream output;
    mmcif_output::MmcifWriter{output}.write_pdb(record, charge_set(1), source);
    auto document = ::gemmi::cif::read_string(output.str());
    auto atom_sites =
        document.blocks.front().find("_atom_site.", {"label_atom_id", "label_alt_id"});
    REQUIRE(atom_sites.length() == 2);
    CHECK(::gemmi::cif::as_string(atom_sites[0][0]) == "C1");
    CHECK(::gemmi::cif::as_string(atom_sites[0][1]) == "A");
    CHECK(::gemmi::cif::as_string(atom_sites[1][0]) == "O1");
}
