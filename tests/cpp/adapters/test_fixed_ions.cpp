#include "adapters/fixed_ions.h"

#include <chargefw/adapters/gemmi/input_options.h>
#include <chargefw/adapters/gemmi/mmcif_input.h>
#include <chargefw/calculation/calculation.h>
#include <chargefw/core/atom.h>
#include <chargefw/core/bond.h>
#include <chargefw/core/molecule.h>
#include <chargefw/core/molecule_collection.h>
#include <chargefw/parameters/io/parameter_set_io.h>
#include <snitch/snitch.hpp>

#include <algorithm>
#include <array>
#include <filesystem>
#include <sstream>
#include <stdexcept>
#include <string>
#include <string_view>
#include <vector>

namespace detail = chargefw::adapters::detail;
using chargefw::adapters::ImportedMoleculeRecord;
namespace calculation = chargefw::calculation;
namespace gemmi_adapter = chargefw::adapters::gemmi;
namespace mmcif = chargefw::adapters::gemmi::mmcif_input;

namespace {
auto mmcif_water_and_magnesium() -> ImportedMoleculeRecord {
    std::istringstream input{R"cif(data_fixed_ions
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
_atom_site.label_entity_id
_atom_site.pdbx_PDB_model_num
HETATM 1 O O . HOH W . ? 0 0 0 1 20 0 5 HOH W O 1 1
HETATM 2 H H1 . HOH W . ? 0.96 0 0 1 20 0 5 HOH W H1 1 1
HETATM 3 H H2 . HOH W . ? -0.24 0.93 0 1 20 0 5 HOH W H2 1 1
HETATM 4 Mg MG . MG I . ? 0 0 3 1 20 0 9 MG I MG 2 1
HETATM 5 O O . HOH W . ? 0 0 0 1 20 0 5 HOH W O 1 2
HETATM 6 H H1 . HOH W . ? 0.96 0 0 1 20 0 5 HOH W H1 1 2
HETATM 7 H H2 . HOH W . ? -0.24 0.93 0 1 20 0 5 HOH W H2 1 2
HETATM 8 Mg MG . MG I . ? 0.25 0 3 1 20 0 9 MG I MG 2 2
#
)cif"};
    auto options = gemmi_adapter::InputOptions{};
    options.bond_strategy = gemmi_adapter::BondStrategy::templates;
    auto reader = mmcif::MmcifReader{input, "fixed-ion-test.cif", options};
    auto record = reader.next();
    if (!record || reader.next()) {
        throw std::runtime_error{"test mmCIF should contain one record"};
    }
    return std::move(*record);
}

auto ion_record(const std::string& component, const int atomic_number, const std::string& atom_name,
                const std::string& sequence, const int formal_charge = 0)
    -> ImportedMoleculeRecord {
    ImportedMoleculeRecord record{.molecule = chargefw::core::Molecule{{chargefw::core::Atom{
                                      atomic_number, formal_charge, atom_name}}}};
    chargefw::adapters::MoleculeImportMetadata metadata;
    metadata.atoms.push_back(
        {.position = 100,
         .structural_labels = chargefw::adapters::SourceStructuralLabels{
             .author =
                 {.atom = atom_name, .residue = component, .chain = " ", .sequence = sequence},
             .label = {
                 .atom = atom_name, .residue = component, .chain = "L", .sequence = sequence}}});
    metadata.components.push_back({component, {0}});
    record.import_metadata = std::move(metadata);
    return record;
}

} // namespace

TEST_CASE("private fixed-charge ion selections share the resolver catalog",
          "[adapters][fixed-charge]") {
    const auto common = detail::fixed_ion_names(true);
    const auto all = detail::fixed_ion_names(false);
    const std::vector<std::string> expected{"NA", "K", "MG", "CA", "CL", "ZN", "FE", "FE2"};
    CHECK(common == expected);
    CHECK(all.size() > common.size());
    for (const auto& name : common) {
        CHECK(std::ranges::count(all, name) == 1);
    }
    for (const auto& name : all) {
        CHECK(std::ranges::count(all, name) == 1);
    }
    const std::vector<ImportedMoleculeRecord> records{ion_record("LIG", 6, "C", "1")};
    CHECK(detail::resolve_fixed_ions(records, all).sources.empty());
}

TEST_CASE("private fixed ions resolve exact ion identities in input order",
          "[adapters][fixed-charge]") {
    const std::vector records{
        ion_record("FE2", 26, "FE", "1", 4), ion_record("MG", 12, "MG", "2", -2),
        ion_record("FE", 26, "FE", "3", 1),  ion_record("NA", 11, "NA", "4", -1),
        ion_record("K", 19, "K", "5", -1),   ion_record("CA", 20, "CA", "6", -1),
        ion_record("CL", 17, "CL", "7", 1),  ion_record("ZN", 30, "ZN", "8", -1)};
    const std::vector<std::string> request{"ZN", "FE", "CA", "K", "NA", "CL", "MG", "FE2"};
    const auto fixed_ions = detail::resolve_fixed_ions(records, request);
    REQUIRE(fixed_ions.sources.size() == 8);
    CHECK(fixed_ions.sources[0].molecule_index == 0);
    CHECK(fixed_ions.sources[0].charge == 2.0);
    CHECK(fixed_ions.sources[1].molecule_index == 1);
    CHECK(fixed_ions.sources[1].charge == 2.0);
    CHECK(fixed_ions.sources[2].molecule_index == 2);
    CHECK(fixed_ions.sources[2].charge == 3.0);
    CHECK(fixed_ions.sources[3].charge == 1.0);
    CHECK(fixed_ions.sources[4].charge == 1.0);
    CHECK(fixed_ions.sources[5].charge == 2.0);
    CHECK(fixed_ions.sources[6].charge == -1.0);
    CHECK(fixed_ions.sources[7].charge == 2.0);
    for (std::size_t i = 0; i < records.size(); ++i) {
        CHECK(fixed_ions.sources[i].atom_index == 0);
    }
    for (std::size_t i = 0; i < records.size(); ++i) {
        CHECK(records[i].molecule.atom(0).formal_charge() != fixed_ions.sources[i].charge);
    }
    CHECK(records[0].molecule.atom(0).formal_charge() == 4);
}

TEST_CASE("private fixed-charge catalog resolves verified CCD monatomic ions",
          "[adapters][fixed-charge]") {
    struct ExpectedIon {
        const char* component;
        const char* atom;
        int element;
        double charge;
    };
    // Independent expected values from the verified 2026-10-05 CCD snapshot.
    const std::array expected{
        ExpectedIon{"0BE", "BE", 4, 2},  ExpectedIon{"3CO", "CO", 27, 3},
        ExpectedIon{"3NI", "NI", 28, 3}, ExpectedIon{"4MO", "MO", 42, 4},
        ExpectedIon{"4PU", "PU", 94, 4}, ExpectedIon{"4TI", "TI", 22, 4},
        ExpectedIon{"6MO", "MO", 42, 6}, ExpectedIon{"AG", "AG", 47, 1},
        ExpectedIon{"AL", "AL", 13, 3},  ExpectedIon{"AM", "AM", 95, 3},
        ExpectedIon{"AU", "AU", 79, 1},  ExpectedIon{"AU3", "AU", 79, 3},
        ExpectedIon{"BA", "BA", 56, 2},  ExpectedIon{"BR", "BR", 35, -1},
        ExpectedIon{"BS3", "BI", 83, 3}, ExpectedIon{"CA", "CA", 20, 2},
        ExpectedIon{"CD", "CD", 48, 2},  ExpectedIon{"CE", "CE", 58, 3},
        ExpectedIon{"CF", "CF", 98, 3},  ExpectedIon{"CL", "CL", 17, -1},
        ExpectedIon{"CO", "CO", 27, 2},  ExpectedIon{"CR", "CR", 24, 3},
        ExpectedIon{"CS", "CS", 55, 1},  ExpectedIon{"CU", "CU", 29, 2},
        ExpectedIon{"CU1", "CU", 29, 1}, ExpectedIon{"CU3", "CU", 29, 3},
        ExpectedIon{"D8U", "D", 1, 1},   ExpectedIon{"DY", "DY", 66, 3},
        ExpectedIon{"ER3", "ER", 68, 3}, ExpectedIon{"EU", "EU", 63, 2},
        ExpectedIon{"EU3", "EU", 63, 3}, ExpectedIon{"F", "F", 9, -1},
        ExpectedIon{"FE", "FE", 26, 3},  ExpectedIon{"FE2", "FE", 26, 2},
        ExpectedIon{"GA", "GA", 31, 3},  ExpectedIon{"GD3", "GD", 64, 3},
        ExpectedIon{"HG", "HG", 80, 2},  ExpectedIon{"HO3", "HO", 67, 3},
        ExpectedIon{"IN", "IN", 49, 3},  ExpectedIon{"IOD", "I", 53, -1},
        ExpectedIon{"IR", "IR", 77, 4},  ExpectedIon{"IR3", "IR", 77, 3},
        ExpectedIon{"K", "K", 19, 1},    ExpectedIon{"LA", "LA", 57, 3},
        ExpectedIon{"LI", "LI", 3, 1},   ExpectedIon{"LU", "LU", 71, 3},
        ExpectedIon{"MG", "MG", 12, 2},  ExpectedIon{"MN", "MN", 25, 2},
        ExpectedIon{"MN3", "MN", 25, 3}, ExpectedIon{"NA", "NA", 11, 1},
        ExpectedIon{"ND", "ND", 60, 3},  ExpectedIon{"NI", "NI", 28, 2},
        ExpectedIon{"OS", "OS", 76, 3},  ExpectedIon{"OS4", "OS", 76, 4},
        ExpectedIon{"PB", "PB", 82, 2},  ExpectedIon{"PD", "PD", 46, 2},
        ExpectedIon{"PR", "PR", 59, 3},  ExpectedIon{"PT", "PT", 78, 2},
        ExpectedIon{"PT4", "PT", 78, 4}, ExpectedIon{"RB", "RB", 37, 1},
        ExpectedIon{"RH", "RH1", 45, 1}, ExpectedIon{"RH3", "RH", 45, 3},
        ExpectedIon{"RHF", "RH", 45, 2}, ExpectedIon{"RU", "RU", 44, 3},
        ExpectedIon{"SB", "SB", 51, 3},  ExpectedIon{"SM", "SM", 62, 3},
        ExpectedIon{"SR", "SR", 38, 2},  ExpectedIon{"TB", "TB", 65, 3},
        ExpectedIon{"TH", "TH", 90, 4},  ExpectedIon{"TL", "TL", 81, 1},
        ExpectedIon{"V", "V", 23, 3},    ExpectedIon{"W", "W", 74, 6},
        ExpectedIon{"Y1", "Y", 39, 2},   ExpectedIon{"YB", "YB", 70, 3},
        ExpectedIon{"YB2", "YB", 70, 2}, ExpectedIon{"YT3", "Y", 39, 3},
        ExpectedIon{"ZCM", "CM", 96, 3}, ExpectedIon{"ZN", "ZN", 30, 2},
        ExpectedIon{"ZR", "ZR", 40, 4},  ExpectedIon{"ZTM", "AC", 89, 3},
    };
    const auto all = detail::fixed_ion_names(false);
    REQUIRE(all.size() == expected.size());
    std::vector<ImportedMoleculeRecord> records;
    auto expected_ion = expected.begin();
    for (const auto& name : all) {
        const auto& ion = *expected_ion++;
        CHECK(name == ion.component);
        records.push_back(ion_record(ion.component, ion.element, ion.atom, "1"));
    }
    for (const auto& ion : {ExpectedIon{"MG", "MG", 12, 2}, ExpectedIon{"FE2", "FE", 26, 2},
                            ExpectedIon{"D8U", "D", 1, 1}}) {
        CAPTURE(ion.component);
        const std::vector<std::string> selected{ion.component};
        const std::vector wrong_name{ion_record(ion.component, ion.element, "WRONG", "1")};
        CHECK(detail::resolve_fixed_ions(wrong_name, selected).sources[0].charge == ion.charge);
        const std::vector wrong_element{
            ion_record(ion.component, ion.element == 1 ? 2 : 1, ion.atom, "1")};
        CHECK_THROWS_AS(detail::resolve_fixed_ions(wrong_element, selected), std::invalid_argument);
        if (std::string_view{ion.atom} != ion.component) {
            const std::vector component_as_atom{
                ion_record(ion.component, ion.element, ion.component, "1")};
            CHECK(detail::resolve_fixed_ions(component_as_atom, selected).sources[0].charge ==
                  ion.charge);
        }
    }
    const auto resolved = detail::resolve_fixed_ions(records, all);
    REQUIRE(resolved.sources.size() == expected.size());
    expected_ion = expected.begin();
    for (std::size_t i = 0; i < expected.size(); ++i) {
        CHECK(resolved.sources[i].molecule_index == i);
        CHECK(resolved.sources[i].atom_index == 0);
        CHECK(resolved.sources[i].charge == expected_ion->charge);
        CHECK(records[i].molecule.atom(0).formal_charge() == 0);
        ++expected_ion;
    }
    for (const std::string excluded : {"ZN2", "SO4", "XE"}) {
        const std::vector<std::string> selected{excluded};
        CHECK_THROWS_AS(detail::resolve_fixed_ions(records, selected), std::invalid_argument);
    }
}

TEST_CASE("private fixed ions resolve repeated instances in atom order",
          "[adapters][fixed-charge]") {
    auto record = ion_record("MG", 12, "MG", "1");
    record.molecule = chargefw::core::Molecule{
        {chargefw::core::Atom{12, 0, "MG"}, chargefw::core::Atom{12, 0, "MG"}}};
    auto second = record.import_metadata->atoms.front();
    second.position = 40;
    record.import_metadata->atoms.push_back(second);
    record.import_metadata->components = {{"MG", {1}}, {"MG", {0}}};

    const std::vector<std::string> request{"MG", "MG"};
    const auto fixed_ions = detail::resolve_fixed_ions(std::vector{record}, request);
    REQUIRE(fixed_ions.sources.size() == 2);
    CHECK(fixed_ions.sources[0].atom_index == 0);
    CHECK(fixed_ions.sources[1].atom_index == 1);
}

TEST_CASE("fixed calcium recognition uses component identity, element, and cardinality",
          "[adapters][fixed-charge]") {
    const std::vector<std::string> selected{"CA"};
    auto calcium = ion_record("CA", 20, "arbitrary-site-name", "1", -1);
    calcium.import_metadata->atoms[0].structural_labels.reset();
    const auto resolved = detail::resolve_fixed_ions(std::vector{calcium}, selected);
    REQUIRE(resolved.sources.size() == 1);
    CHECK(resolved.sources[0].atom_index == 0);
    CHECK(resolved.sources[0].charge == 2.0);
    CHECK(calcium.molecule.atom(0).formal_charge() == -1);

    const auto singleton = ion_record("LIG", 20, "CA", "1");
    CHECK(detail::resolve_fixed_ions(std::vector{singleton}, selected).sources.empty());
    const auto carbon = ion_record("CA", 6, "CA", "1");
    CHECK_THROWS_AS(detail::resolve_fixed_ions(std::vector{carbon}, selected),
                    std::invalid_argument);

    calcium.molecule =
        chargefw::core::Molecule{{chargefw::core::Atom{20}, chargefw::core::Atom{20}}};
    calcium.import_metadata->atoms.push_back(calcium.import_metadata->atoms.front());
    calcium.import_metadata->components[0].atom_indices.push_back(1);
    CHECK_THROWS_AS(detail::resolve_fixed_ions(std::vector{calcium}, selected),
                    std::invalid_argument);
}

TEST_CASE("private fixed ions preserve atom order within one record", "[adapters][fixed-charge]") {
    auto record = ion_record("FE2", 26, "FE", "9");
    record.molecule = chargefw::core::Molecule{{chargefw::core::Atom{26, 0, "FE"},
                                                chargefw::core::Atom{12, 7, "MG"},
                                                chargefw::core::Atom{26, -4, "FE"}}};
    auto first = record.import_metadata->atoms.front();
    auto second = first;
    auto third = first;
    first.position = 900;
    second.position = 2;
    third.position = 40;
    first.structural_labels->label.sequence = "9";
    first.structural_labels->label.atom = "FE";
    first.structural_labels->author.atom = "FE";
    second.structural_labels->label.sequence = "1";
    second.structural_labels->author.sequence = "1";
    second.structural_labels->label.residue = "MG";
    second.structural_labels->author.residue = "MG";
    second.structural_labels->label.atom = "MG";
    second.structural_labels->author.atom = "MG";
    third.structural_labels->label.sequence = "5";
    third.structural_labels->author.sequence = "5";
    third.structural_labels->label.residue = "FE";
    third.structural_labels->author.residue = "FE";
    third.structural_labels->label.atom = "FE";
    third.structural_labels->author.atom = "FE";
    record.import_metadata->atoms = {first, second, third};
    record.import_metadata->components = {{"FE2", {0}}, {"MG", {1}}, {"FE", {2}}};

    const std::vector records{record};
    const std::vector<std::string> request{"FE", "MG", "FE2", "MG"};
    const auto fixed_ions = detail::resolve_fixed_ions(records, request);
    REQUIRE(fixed_ions.sources.size() == 3);
    CHECK(fixed_ions.sources[0].atom_index == 0);
    CHECK(fixed_ions.sources[0].charge == 2.0);
    CHECK(fixed_ions.sources[1].atom_index == 1);
    CHECK(fixed_ions.sources[1].charge == 2.0);
    CHECK(fixed_ions.sources[2].atom_index == 2);
    CHECK(fixed_ions.sources[2].charge == 3.0);
    CHECK(records[0].molecule.atom(0).formal_charge() == 0);
    CHECK(records[0].molecule.atom(1).formal_charge() == 7);
    CHECK(records[0].molecule.atom(2).formal_charge() == -4);
}

TEST_CASE("private fixed ions handle empty, absent, and unknown selections",
          "[adapters][fixed-charge]") {
    const std::vector<ImportedMoleculeRecord> records;
    const std::vector<std::string> none;
    CHECK(detail::resolve_fixed_ions(records, none).sources.empty());
    const auto no_metadata =
        ImportedMoleculeRecord{.molecule = chargefw::core::Molecule{{chargefw::core::Atom{6}}}};
    const std::vector malformed{no_metadata};
    CHECK(detail::resolve_fixed_ions(malformed, none).sources.empty());
    const std::vector<std::string> absent{"NA"};
    const std::vector valid{ion_record("LIG", 6, "C", "1")};
    const auto absent_result = detail::resolve_fixed_ions(valid, absent);
    CHECK(absent_result.sources.empty());
    for (const std::string bad : {"", "mg", "Mg"}) {
        const std::vector<std::string> unknown{bad};
        CHECK_THROWS_AS(detail::resolve_fixed_ions(records, unknown), std::invalid_argument);
    }
}

TEST_CASE("private fixed ions use explicit components without hierarchy labels",
          "[adapters][fixed-charge]") {
    auto record = ion_record("MG", 12, "MG", "1");
    record.import_metadata->atoms[0].structural_labels.reset();
    const std::vector records{record};
    const std::vector<std::string> request{"MG"};
    const auto fixed_ions = detail::resolve_fixed_ions(records, request);
    REQUIRE(fixed_ions.sources.size() == 1);
    CHECK(fixed_ions.sources[0].charge == 2.0);
}

TEST_CASE("private fixed ions use canonical component IDs and validate partitions",
          "[adapters][fixed-charge]") {
    auto record = ion_record("MG", 12, "MG", "1");
    record.import_metadata->atoms[0].structural_labels->author.residue = "AUTHOR-NAMESPACE";
    record.import_metadata->components = {{"MG", {0}}};
    const std::vector<std::string> request{"MG"};
    const auto resolved = detail::resolve_fixed_ions(std::vector{record}, request);
    REQUIRE(resolved.sources.size() == 1);
    CHECK(resolved.sources[0].atom_index == 0);

    auto malformed = record;
    malformed.import_metadata->components = {{"MG", {1}}};
    CHECK_THROWS_AS(detail::resolve_fixed_ions(std::vector{malformed}, request),
                    std::invalid_argument);
    malformed.import_metadata->components = {{"MG", {0}}, {"MG", {0}}};
    CHECK_THROWS_AS(detail::resolve_fixed_ions(std::vector{malformed}, request),
                    std::invalid_argument);
    malformed.import_metadata->components = {{"", {0}}};
    CHECK_THROWS_AS(detail::resolve_fixed_ions(std::vector{malformed}, request),
                    std::invalid_argument);

    auto incomplete = record;
    incomplete.molecule = chargefw::core::Molecule{
        {chargefw::core::Atom{12, 0, "MG"}, chargefw::core::Atom{1, 0, "H"}}};
    auto second = incomplete.import_metadata->atoms.front();
    second.position = 2;
    incomplete.import_metadata->atoms.push_back(second);
    incomplete.import_metadata->components = {{"MG", {0}}};
    CHECK_THROWS_AS(detail::resolve_fixed_ions(std::vector{incomplete}, request),
                    std::invalid_argument);
}

TEST_CASE("private fixed ions match label names before author fallback",
          "[adapters][fixed-charge]") {
    auto record = ion_record("MG", 12, "AUTHOR-MG", "1");
    record.import_metadata->atoms[0].structural_labels->author.atom = "AUTHOR-MG";
    record.import_metadata->atoms[0].structural_labels->label.atom = "MG";
    record.import_metadata->atoms[0].structural_labels->author.residue = "AUTHOR-COMP";
    record.import_metadata->atoms[0].structural_labels->label.residue = "MG";
    const std::vector<std::string> request{"MG"};
    const auto resolved = detail::resolve_fixed_ions(std::vector{record}, request);
    REQUIRE(resolved.sources.size() == 1);
    CHECK(resolved.sources[0].atom_index == 0);

    record.import_metadata->atoms[0].structural_labels->label.residue = "NOT-MG";
    CHECK(detail::resolve_fixed_ions(std::vector{record}, request).sources.empty());
    record.import_metadata->atoms[0].structural_labels->label.residue = "";
    record.import_metadata->atoms[0].structural_labels->author.residue = "MG";
    CHECK(detail::resolve_fixed_ions(std::vector{record}, request).sources.size() == 1);
    record.import_metadata->components.clear();
    CHECK_THROWS_AS(detail::resolve_fixed_ions(std::vector{record}, request),
                    std::invalid_argument);
}

TEST_CASE("private fixed ions reject incomplete and inconsistent mappings",
          "[adapters][fixed-charge]") {
    const std::vector<std::string> request{"MG"};
    const auto missing =
        ImportedMoleculeRecord{.molecule = chargefw::core::Molecule{{chargefw::core::Atom{12}}}};
    CHECK_THROWS_AS(detail::resolve_fixed_ions(std::vector{missing}, request),
                    std::invalid_argument);

    auto wrong_size = ion_record("MG", 12, "MG", "1");
    wrong_size.import_metadata->atoms.clear();
    CHECK_THROWS_AS(detail::resolve_fixed_ions(std::vector{wrong_size}, request),
                    std::invalid_argument);
}

TEST_CASE("private fixed ions reject extra atoms and incident bonds", "[adapters][fixed-charge]") {
    const std::vector<std::string> mg_request{"MG"};
    auto extra_atom = ion_record("MG", 12, "MG", "1");
    extra_atom.molecule = chargefw::core::Molecule{
        {chargefw::core::Atom{12, 0, "MG"}, chargefw::core::Atom{1, 0, "H"}}};
    auto extra_metadata = extra_atom.import_metadata->atoms.front();
    extra_metadata.position = 101;
    extra_metadata.structural_labels->author.atom = "H";
    extra_metadata.structural_labels->label.atom = "H";
    extra_atom.import_metadata->atoms.push_back(extra_metadata);
    extra_atom.import_metadata->components[0].atom_indices.push_back(1);
    CHECK_THROWS_AS(detail::resolve_fixed_ions(std::vector{extra_atom}, mg_request),
                    std::invalid_argument);

    auto bonded = ion_record("MG", 12, "MG", "1");
    bonded.molecule = chargefw::core::Molecule{
        {chargefw::core::Atom{12, 0, "MG"}, chargefw::core::Atom{8, 0, "O"}},
        {chargefw::core::Bond{0, 1}}};
    auto active_metadata = bonded.import_metadata->atoms.front();
    active_metadata.position = 101;
    active_metadata.structural_labels->author.atom = "O";
    active_metadata.structural_labels->author.residue = "HOH";
    active_metadata.structural_labels->author.sequence = "2";
    active_metadata.structural_labels->label.atom = "O";
    active_metadata.structural_labels->label.residue = "HOH";
    active_metadata.structural_labels->label.sequence = "2";
    bonded.import_metadata->atoms.push_back(active_metadata);
    bonded.import_metadata->components.push_back({"HOH", {1}});
    CHECK_THROWS_AS(detail::resolve_fixed_ions(std::vector{bonded}, mg_request),
                    std::invalid_argument);

    auto two_selected = ion_record("MG", 12, "MG", "1");
    two_selected.molecule = chargefw::core::Molecule{
        {chargefw::core::Atom{12, 0, "MG"}, chargefw::core::Atom{12, 0, "MG"}},
        {chargefw::core::Bond{0, 1}}};
    auto second = two_selected.import_metadata->atoms.front();
    second.position = 101;
    second.structural_labels->author.sequence = "2";
    second.structural_labels->label.sequence = "2";
    two_selected.import_metadata->atoms.push_back(second);
    two_selected.import_metadata->components = {{"MG", {0}}, {"MG", {1}}};
    auto bond_error = std::string{};
    try {
        static_cast<void>(detail::resolve_fixed_ions(std::vector{two_selected}, mg_request));
    } catch (const std::invalid_argument& error) {
        bond_error = error.what();
    }
    CHECK(bond_error.contains("incident graph bond"));
}

TEST_CASE("mmCIF named-ion resolution feeds SQE+qp full calculation",
          "[adapters][fixed-charge][sqeqp][mmcif]") {
    const auto record = mmcif_water_and_magnesium();
    REQUIRE(record.molecule.atom_count() == 4);
    REQUIRE(record.molecule.bond_count() == 2);
    REQUIRE(record.molecule.conformer_count() == 2);
    REQUIRE(record.import_metadata.has_value());
    REQUIRE(record.import_metadata->components.size() == 2);
    CHECK(record.import_metadata->components[0] ==
          chargefw::adapters::SourceComponentInstance{"HOH", {0, 1, 2}});
    CHECK(record.import_metadata->components[1] ==
          chargefw::adapters::SourceComponentInstance{"MG", {3}});
    CHECK(record.import_metadata->atoms[3].structural_labels->author.sequence == "9");

    const std::vector records{record};
    const std::vector<std::string> names{"MG"};
    const auto resolved = detail::resolve_fixed_ions(records, names);
    REQUIRE(resolved.sources.size() == 1);
    CHECK(resolved.sources[0].molecule_index == 0);
    CHECK(resolved.sources[0].atom_index == 3);
    CHECK(resolved.sources[0].charge == 2.0);

    const auto original = chargefw::core::MoleculeCollection{{record.molecule}};
    const auto parameter_set = chargefw::parameters::load_parameter_set_json_file(
        std::filesystem::path{CHARGEFW_TEST_PARAMETER_DIR} / "SQEqp_Schindler2021_CCD_gen.json");

    const auto unsupported_assessment =
        calculation::assess({.molecules = original,
                             .parameter_sets = {parameter_set},
                             .method_id = "sqeqp",
                             .parameter_set_id = "SQEqp_Schindler2021_CCD_gen"});
    CHECK(unsupported_assessment.plans().empty());

    const auto fixed_source_result = calculation::calculate(calculation::assess(
        {.molecules = original,
         .parameter_sets = {parameter_set},
         .method_id = "sqeqp",
         .parameter_set_id = "SQEqp_Schindler2021_CCD_gen",
         .execution_selection =
             calculation::ExecutionSelection{calculation::ExecutionSelectionKind::full},
         .fixed_ions = resolved}));
    REQUIRE(fixed_source_result.calculated());
    REQUIRE(fixed_source_result.charges->size() == 2);
    for (std::size_t conformer = 0; conformer < 2; ++conformer) {
        const auto& got = fixed_source_result.charges->assignment(conformer);
        CHECK(got.target.conformer_index == conformer);
        REQUIRE(got.charges.size() == 4);
        CHECK(got.charges[3] == 2.0);
    }
    CHECK(record.molecule.atom(3).formal_charge() == 0);
    CHECK(record.molecule.conformer(0).positions()[3].x == 0.0);
    CHECK(record.molecule.conformer(1).positions()[3].x == 0.25);
}
