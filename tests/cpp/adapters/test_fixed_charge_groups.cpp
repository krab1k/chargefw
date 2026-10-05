#include "adapters/fixed_charge_groups.h"

#include <chargefw/adapters/gemmi/input_options.h>
#include <chargefw/adapters/gemmi/mmcif_input.h>
#include <chargefw/calculation/calculation.h>
#include <chargefw/core/atom.h>
#include <chargefw/core/bond.h>
#include <chargefw/core/conformer.h>
#include <chargefw/core/molecule.h>
#include <chargefw/core/molecule_collection.h>
#include <chargefw/parameters/io/parameter_set_io.h>
#include <snitch/snitch.hpp>

#include <array>
#include <cmath>
#include <filesystem>
#include <optional>
#include <sstream>
#include <stdexcept>
#include <string>
#include <vector>

namespace detail = chargefw::adapters::detail;
using chargefw::adapters::ImportedMoleculeRecord;
namespace calculation = chargefw::calculation;
namespace gemmi_adapter = chargefw::adapters::gemmi;
namespace mmcif = chargefw::adapters::gemmi::mmcif_input;

namespace {
auto mmcif_water_and_magnesium() -> ImportedMoleculeRecord {
    std::istringstream input{R"cif(data_embedding
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

auto run_sqeqp(const chargefw::core::MoleculeCollection& molecules,
               const chargefw::parameters::ParameterSet& parameter_set,
               const calculation::ExecutionSelection execution,
               const std::optional<calculation::FixedChargeEmbedding>& embedding = {})
    -> calculation::ExecutionResult {
    auto request = calculation::AssessmentRequest{.molecules = molecules,
                                                  .parameter_sets = {parameter_set},
                                                  .method_id = "sqeqp",
                                                  .parameter_set_id = "SQEqp_Schindler2021_CCD_gen",
                                                  .execution_selection = execution,
                                                  .fixed_charge_embedding = embedding};
    const auto assessment = calculation::assess(std::move(request));
    return calculation::calculate(assessment);
}
} // namespace

TEST_CASE("private fixed-charge groups resolve exact ion identities in input order",
          "[adapters][fixed-charge]") {
    const std::vector records{
        ion_record("FE2", 26, "FE2", "1", 4), ion_record("MG", 12, "MG", "2", -2),
        ion_record("FE", 26, "FE", "3", 1),   ion_record("NA", 11, "NA", "4", -1),
        ion_record("K", 19, "K", "5", -1),    ion_record("CA", 20, "CA", "6", -1),
        ion_record("CL", 17, "CL", "7", 1),   ion_record("ZN", 30, "ZN", "8", -1)};
    const std::vector<std::string> request{"ZN", "FE", "CA", "K", "NA", "CL", "MG", "FE2"};
    const auto embedding = detail::resolve_fixed_charge_groups(records, request);
    REQUIRE(embedding.sources.size() == 8);
    CHECK(embedding.sources[0].molecule_index == 0);
    CHECK(embedding.sources[0].charge == 2.0);
    CHECK(embedding.sources[1].molecule_index == 1);
    CHECK(embedding.sources[1].charge == 2.0);
    CHECK(embedding.sources[2].molecule_index == 2);
    CHECK(embedding.sources[2].charge == 3.0);
    CHECK(embedding.sources[3].charge == 1.0);
    CHECK(embedding.sources[4].charge == 1.0);
    CHECK(embedding.sources[5].charge == 2.0);
    CHECK(embedding.sources[6].charge == -1.0);
    CHECK(embedding.sources[7].charge == 2.0);
    for (std::size_t i = 0; i < records.size(); ++i) {
        CHECK(embedding.sources[i].atom_index == 0);
    }
    CHECK(embedding.charge_provenance == "chargefw:fixed-charge-ions:v1");
    for (std::size_t i = 0; i < records.size(); ++i) {
        CHECK(records[i].molecule.atom(0).formal_charge() != embedding.sources[i].charge);
    }
    CHECK(records[0].molecule.atom(0).formal_charge() == 4);
}

TEST_CASE("private fixed-charge groups resolve repeated instances in atom order",
          "[adapters][fixed-charge]") {
    auto record = ion_record("MG", 12, "MG", "1");
    record.molecule = chargefw::core::Molecule{
        {chargefw::core::Atom{12, 0, "MG"}, chargefw::core::Atom{12, 0, "MG"}}};
    auto second = record.import_metadata->atoms.front();
    second.position = 40;
    record.import_metadata->atoms.push_back(second);
    record.import_metadata->components = {{"MG", {1}}, {"MG", {0}}};

    const std::vector<std::string> request{"MG", "MG"};
    const auto embedding = detail::resolve_fixed_charge_groups(std::vector{record}, request);
    REQUIRE(embedding.sources.size() == 2);
    CHECK(embedding.sources[0].atom_index == 0);
    CHECK(embedding.sources[1].atom_index == 1);
}

TEST_CASE("private fixed-charge groups preserve atom order within one record",
          "[adapters][fixed-charge]") {
    auto record = ion_record("FE2", 26, "FE2", "9");
    record.molecule = chargefw::core::Molecule{{chargefw::core::Atom{26, 0, "FE2"},
                                                chargefw::core::Atom{12, 7, "MG"},
                                                chargefw::core::Atom{26, -4, "FE"}}};
    auto first = record.import_metadata->atoms.front();
    auto second = first;
    auto third = first;
    first.position = 900;
    second.position = 2;
    third.position = 40;
    first.structural_labels->label.sequence = "9";
    first.structural_labels->label.atom = "FE2";
    first.structural_labels->author.atom = "FE2";
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
    const auto embedding = detail::resolve_fixed_charge_groups(records, request);
    REQUIRE(embedding.sources.size() == 3);
    CHECK(embedding.sources[0].atom_index == 0);
    CHECK(embedding.sources[0].charge == 2.0);
    CHECK(embedding.sources[1].atom_index == 1);
    CHECK(embedding.sources[1].charge == 2.0);
    CHECK(embedding.sources[2].atom_index == 2);
    CHECK(embedding.sources[2].charge == 3.0);
    CHECK(records[0].molecule.atom(0).formal_charge() == 0);
    CHECK(records[0].molecule.atom(1).formal_charge() == 7);
    CHECK(records[0].molecule.atom(2).formal_charge() == -4);
}

TEST_CASE("private fixed-charge groups handle empty, absent, and unknown selections",
          "[adapters][fixed-charge]") {
    const std::vector<ImportedMoleculeRecord> records;
    const std::vector<std::string> none;
    CHECK(detail::resolve_fixed_charge_groups(records, none).sources.empty());
    const auto no_metadata =
        ImportedMoleculeRecord{.molecule = chargefw::core::Molecule{{chargefw::core::Atom{6}}}};
    const std::vector malformed{no_metadata};
    CHECK(detail::resolve_fixed_charge_groups(malformed, none).sources.empty());
    const std::vector<std::string> absent{"NA"};
    const std::vector valid{ion_record("LIG", 6, "C", "1")};
    const auto absent_result = detail::resolve_fixed_charge_groups(valid, absent);
    CHECK(absent_result.sources.empty());
    CHECK(absent_result.charge_provenance.empty());
    for (const std::string bad : {"", "mg", "Mg"}) {
        const std::vector<std::string> unknown{bad};
        CHECK_THROWS_AS(detail::resolve_fixed_charge_groups(records, unknown),
                        std::invalid_argument);
    }
}

TEST_CASE("private fixed-charge groups require isolated matching atom identities",
          "[adapters][fixed-charge]") {
    const std::vector<std::string> request{"MG"};
    auto malformed = ion_record("MG", 12, "X", "1");
    const std::vector records{malformed};
    CHECK_THROWS_AS(detail::resolve_fixed_charge_groups(records, request), std::invalid_argument);
}

TEST_CASE("private fixed-charge groups use explicit components without hierarchy labels",
          "[adapters][fixed-charge]") {
    auto record = ion_record("MG", 12, "MG", "1");
    record.import_metadata->atoms[0].structural_labels.reset();
    const std::vector records{record};
    const std::vector<std::string> request{"MG"};
    const auto embedding = detail::resolve_fixed_charge_groups(records, request);
    REQUIRE(embedding.sources.size() == 1);
    CHECK(embedding.sources[0].charge == 2.0);
}

TEST_CASE("private fixed-charge groups use canonical component IDs and validate partitions",
          "[adapters][fixed-charge]") {
    auto record = ion_record("MG", 12, "MG", "1");
    record.import_metadata->atoms[0].structural_labels->author.residue = "AUTHOR-NAMESPACE";
    record.import_metadata->components = {{"MG", {0}}};
    const std::vector<std::string> request{"MG"};
    const auto resolved = detail::resolve_fixed_charge_groups(std::vector{record}, request);
    REQUIRE(resolved.sources.size() == 1);
    CHECK(resolved.sources[0].atom_index == 0);

    auto malformed = record;
    malformed.import_metadata->components = {{"MG", {1}}};
    CHECK_THROWS_AS(detail::resolve_fixed_charge_groups(std::vector{malformed}, request),
                    std::invalid_argument);
    malformed.import_metadata->components = {{"MG", {0}}, {"MG", {0}}};
    CHECK_THROWS_AS(detail::resolve_fixed_charge_groups(std::vector{malformed}, request),
                    std::invalid_argument);
    malformed.import_metadata->components = {{"", {0}}};
    CHECK_THROWS_AS(detail::resolve_fixed_charge_groups(std::vector{malformed}, request),
                    std::invalid_argument);

    auto incomplete = record;
    incomplete.molecule = chargefw::core::Molecule{
        {chargefw::core::Atom{12, 0, "MG"}, chargefw::core::Atom{1, 0, "H"}}};
    auto second = incomplete.import_metadata->atoms.front();
    second.position = 2;
    incomplete.import_metadata->atoms.push_back(second);
    incomplete.import_metadata->components = {{"MG", {0}}};
    CHECK_THROWS_AS(detail::resolve_fixed_charge_groups(std::vector{incomplete}, request),
                    std::invalid_argument);
}

TEST_CASE("private fixed-charge groups match label names before author fallback",
          "[adapters][fixed-charge]") {
    auto record = ion_record("MG", 12, "AUTHOR-MG", "1");
    record.import_metadata->atoms[0].structural_labels->author.atom = "AUTHOR-MG";
    record.import_metadata->atoms[0].structural_labels->label.atom = "MG";
    record.import_metadata->atoms[0].structural_labels->author.residue = "AUTHOR-COMP";
    record.import_metadata->atoms[0].structural_labels->label.residue = "MG";
    const std::vector<std::string> request{"MG"};
    const auto resolved = detail::resolve_fixed_charge_groups(std::vector{record}, request);
    REQUIRE(resolved.sources.size() == 1);
    CHECK(resolved.sources[0].atom_index == 0);

    record.import_metadata->atoms[0].structural_labels->label.atom = "NOT-MG";
    CHECK_THROWS_AS(detail::resolve_fixed_charge_groups(std::vector{record}, request),
                    std::invalid_argument);
}

TEST_CASE("private fixed-charge groups reject incomplete and inconsistent mappings",
          "[adapters][fixed-charge]") {
    const std::vector<std::string> request{"MG"};
    const auto missing =
        ImportedMoleculeRecord{.molecule = chargefw::core::Molecule{{chargefw::core::Atom{12}}}};
    CHECK_THROWS_AS(detail::resolve_fixed_charge_groups(std::vector{missing}, request),
                    std::invalid_argument);

    auto wrong_size = ion_record("MG", 12, "MG", "1");
    wrong_size.import_metadata->atoms.clear();
    CHECK_THROWS_AS(detail::resolve_fixed_charge_groups(std::vector{wrong_size}, request),
                    std::invalid_argument);

    auto wrong_element = ion_record("MG", 11, "MG", "1");
    CHECK_THROWS_AS(detail::resolve_fixed_charge_groups(std::vector{wrong_element}, request),
                    std::invalid_argument);
}

TEST_CASE("private fixed-charge groups reject extra atoms and incident bonds",
          "[adapters][fixed-charge]") {
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
    CHECK_THROWS_AS(detail::resolve_fixed_charge_groups(std::vector{extra_atom}, mg_request),
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
    CHECK_THROWS_AS(detail::resolve_fixed_charge_groups(std::vector{bonded}, mg_request),
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
        static_cast<void>(
            detail::resolve_fixed_charge_groups(std::vector{two_selected}, mg_request));
    } catch (const std::invalid_argument& error) {
        bond_error = error.what();
    }
    CHECK(bond_error.contains("incident graph bond"));
}

TEST_CASE("mmCIF named-ion resolution feeds SQE+qp full and whole-radius calculations",
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
    const auto resolved = detail::resolve_fixed_charge_groups(records, names);
    REQUIRE(resolved.sources.size() == 1);
    CHECK(resolved.sources[0].molecule_index == 0);
    CHECK(resolved.sources[0].atom_index == 3);
    CHECK(resolved.sources[0].charge == 2.0);

    auto active_atoms = std::vector<chargefw::core::Atom>{};
    auto active_bonds = std::vector<chargefw::core::Bond>{};
    for (std::size_t i = 0; i < 3; ++i) {
        active_atoms.push_back(record.molecule.atom(i));
    }
    for (const auto& bond : record.molecule.bonds()) {
        active_bonds.push_back(bond);
    }
    auto active_conformers = std::vector<chargefw::core::Conformer>{};
    for (const auto& conformer : record.molecule.conformers()) {
        const auto positions = conformer.positions();
        active_conformers.emplace_back(
            std::vector<chargefw::core::Position>{positions[0], positions[1], positions[2]},
            std::string{conformer.name()});
    }
    const auto active = chargefw::core::MoleculeCollection{
        {chargefw::core::Molecule{active_atoms, active_bonds, active_conformers, "active-water"}}};
    const auto original = chargefw::core::MoleculeCollection{{record.molecule}};
    const auto parameter_set = chargefw::parameters::load_parameter_set_json_file(
        std::filesystem::path{CHARGEFW_TEST_PARAMETER_DIR} / "SQEqp_Schindler2021_CCD_gen.json");

    const auto unembedded =
        calculation::assess({.molecules = original,
                             .parameter_sets = {parameter_set},
                             .method_id = "sqeqp",
                             .parameter_set_id = "SQEqp_Schindler2021_CCD_gen"});
    CHECK(unembedded.plans().empty());

    const auto embedded = run_sqeqp(
        original, parameter_set,
        calculation::ExecutionSelection{calculation::ExecutionSelectionKind::full}, resolved);
    const auto manual = run_sqeqp(
        original, parameter_set,
        calculation::ExecutionSelection{calculation::ExecutionSelectionKind::full},
        calculation::FixedChargeEmbedding{{{0, 3, 2.0}}, "chargefw:fixed-charge-ions:v1"});
    const auto omitted =
        run_sqeqp(active, parameter_set,
                  calculation::ExecutionSelection{calculation::ExecutionSelectionKind::full});
    REQUIRE(embedded.calculated());
    REQUIRE(manual.calculated());
    REQUIRE(omitted.calculated());
    REQUIRE(embedded.charges->size() == 2);
    for (std::size_t conformer = 0; conformer < 2; ++conformer) {
        const auto& got = embedded.charges->assignment(conformer);
        const auto& expected = manual.charges->assignment(conformer);
        CHECK(got.target.conformer_index == conformer);
        REQUIRE(got.charges.size() == 4);
        CHECK(got.charges[3] == 2.0);
        for (std::size_t atom = 0; atom < got.charges.size(); ++atom) {
            CHECK(std::abs(got.charges[atom] - expected.charges[atom]) < 1.0e-12);
        }
        CHECK(std::abs(got.charges[0] - omitted.charges->assignment(conformer).charges[0]) >
              1.0e-6);
    }

    for (const auto kind : {calculation::ExecutionSelectionKind::cutoff,
                            calculation::ExecutionSelectionKind::cover}) {
        const auto reduced = run_sqeqp(original, parameter_set,
                                       calculation::ExecutionSelection{kind, 8.0}, resolved);
        REQUIRE(reduced.calculated());
        REQUIRE(reduced.charges->size() == embedded.charges->size());
        for (std::size_t assignment = 0; assignment < reduced.charges->size(); ++assignment) {
            const auto& reduced_values = reduced.charges->assignment(assignment).charges;
            const auto& full_values = embedded.charges->assignment(assignment).charges;
            REQUIRE(reduced_values.size() == full_values.size());
            for (std::size_t atom = 0; atom < full_values.size(); ++atom) {
                CHECK(std::abs(reduced_values[atom] - full_values[atom]) < 1.0e-10);
            }
            CHECK(reduced_values[3] == 2.0);
        }
    }
    CHECK(record.molecule.atom(3).formal_charge() == 0);
    CHECK(record.molecule.conformer(0).positions()[3].x == 0.0);
    CHECK(record.molecule.conformer(1).positions()[3].x == 0.25);
}
