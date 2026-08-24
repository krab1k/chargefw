#include "support/test_molecules.h"

#include <chargefw/features/prepared_molecule.h>
#include <chargefw/methods/calculation_input.h>
#include <chargefw/methods/method_options.h>
#include <chargefw/methods/method_registry.h>

#include <cmath>
#include <snitch/snitch.hpp>
#include <string>
#include <string_view>
#include <variant>
#include <vector>

namespace features = chargefw::features;
namespace methods = chargefw::methods;

namespace {

auto calculate(const methods::Method& method, const chargefw::core::Molecule& molecule)
    -> chargefw::charges::AtomicCharges {
    const features::PreparedMolecule prepared_molecule{molecule};
    const auto method_options = methods::make_default_options(method.option_schema());

    const methods::CalculationInput input{prepared_molecule, method_options,
                                          chargefw::core::total_formal_charge(molecule)};

    return method.calculate(input);
}

} // namespace

TEST_CASE("built-in method registry exposes all methods with correct metadata",
          "[methods][builtin-methods]") {
    const auto& registry = methods::method_registry();

    const auto method_names = registry.names();

    CHECK((method_names == std::vector<std::string>{"abeem",  "charge2", "delre", "denr",   "dummy",
                                                    "eem",    "eqeq",    "eqeqc", "formal", "gdac",
                                                    "kcm",    "mgc",     "mpeoe", "peoe",   "qeq",
                                                    "sfkeem", "smpqeq",  "sqe",   "sqeq0",  "sqeqp",
                                                    "tsef",   "veem"}));

    for (const auto& name : method_names) {
        CHECK(registry.find(name) != nullptr);
    }

    const auto* dummy = registry.find("dummy");
    const auto* formal = registry.find("formal");
    const auto* veem = registry.find("veem");
    const auto* peoe = registry.find("peoe");
    const auto* mpeoe = registry.find("mpeoe");
    const auto* gdac = registry.find("gdac");
    const auto* charge2 = registry.find("charge2");
    const auto* delre = registry.find("delre");
    const auto* mgc = registry.find("mgc");
    const auto* denr = registry.find("denr");
    const auto* kcm = registry.find("kcm");
    const auto* tsef = registry.find("tsef");
    const auto* qeq = registry.find("qeq");
    const auto* sqe = registry.find("sqe");
    const auto* sqeq0 = registry.find("sqeq0");
    const auto* sqeqp = registry.find("sqeqp");
    const auto* eem = registry.find("eem");
    const auto* smpqeq = registry.find("smpqeq");
    const auto* sfkeem = registry.find("sfkeem");
    const auto* eqeq = registry.find("eqeq");
    const auto* eqeqc = registry.find("eqeqc");
    const auto* abeem = registry.find("abeem");

    CHECK(dummy->id() == std::string_view{"dummy"});
    CHECK(dummy->metadata().name == std::string_view{"Dummy method"});
    CHECK(dummy->metadata().full_name == std::string_view{"Dummy zero charges"});
    CHECK(!dummy->metadata().publication.has_value());
    CHECK(dummy->metadata().priority == 0);
    CHECK(!dummy->requires_parameters());
    CHECK(dummy->option_schema().empty());

    CHECK(formal->id() == std::string_view{"formal"});
    CHECK(formal->metadata().name == std::string_view{"Formal"});
    CHECK(formal->metadata().full_name == std::string_view{"Formal atomic charges"});
    CHECK(!formal->metadata().publication.has_value());
    CHECK(formal->metadata().priority == 10);
    CHECK(!formal->requires_parameters());
    CHECK(formal->option_schema().empty());

    CHECK(veem->id() == std::string_view{"veem"});
    CHECK(veem->metadata().name == std::string_view{"VEEM"});
    CHECK(veem->metadata().full_name == std::string_view{"Valence Electrons Equalization Method"});
    CHECK(veem->metadata().publication.has_value());
    CHECK(veem->metadata().priority == 20);
    CHECK(veem->requirements().resources.time == methods::ComplexityTerm::atoms);
    CHECK(veem->requirements().resources.memory == methods::ComplexityTerm::constant);
    CHECK(!veem->requires_parameters());
    CHECK(veem->option_schema().empty());

    CHECK(peoe->id() == std::string_view{"peoe"});
    CHECK(peoe->metadata().name == std::string_view{"PEOE"});
    CHECK(peoe->metadata().full_name ==
          std::string_view{"Partial Equalization of Atomic Electronegativity"});
    CHECK(peoe->metadata().publication.has_value());
    CHECK(peoe->metadata().priority == 120);

    CHECK(mpeoe->id() == std::string_view{"mpeoe"});
    CHECK(mpeoe->metadata().name == std::string_view{"MPEOE"});
    CHECK(mpeoe->metadata().full_name ==
          std::string_view{"Modified Partial Equalization of Atomic Electronegativity"});
    CHECK(mpeoe->metadata().publication.has_value());
    CHECK(mpeoe->metadata().priority == 110);

    CHECK(gdac->id() == std::string_view{"gdac"});
    CHECK(gdac->metadata().name == std::string_view{"GDAC"});
    CHECK(gdac->metadata().full_name == std::string_view{"Geometry-Dependent Net Atomic Charges"});
    CHECK(gdac->metadata().publication.has_value());
    CHECK(gdac->metadata().priority == 100);

    CHECK(gdac->requirements().coordinates);
    CHECK(gdac->requirements().requires_atom_parameters());
    CHECK(!gdac->requirements().requires_common_parameters());
    CHECK(!gdac->requirements().requires_bond_parameters());
    CHECK(gdac->requirements().atom_parameters.size() == 2);
    CHECK(gdac->requires_parameters());
    CHECK(gdac->option_schema().size() == 1);

    CHECK(mpeoe->requirements().requires_common_parameters());
    CHECK(mpeoe->requirements().requires_atom_parameters());
    CHECK(mpeoe->requirements().requires_bond_parameters());
    CHECK(mpeoe->requirements().common_parameters.size() == 1);
    CHECK(mpeoe->requirements().atom_parameters.size() == 2);
    CHECK(mpeoe->requirements().bond_parameters.size() == 1);
    CHECK(mpeoe->requires_parameters());
    CHECK(mpeoe->option_schema().size() == 1);

    CHECK(peoe->requirements().requires_common_parameters());
    CHECK(peoe->requirements().requires_atom_parameters());
    CHECK(peoe->requirements().common_parameters.size() == 1);
    CHECK(peoe->requirements().atom_parameters.size() == 3);
    CHECK(peoe->requires_parameters());
    CHECK(peoe->option_schema().size() == 1);

    CHECK(charge2->id() == std::string_view{"charge2"});
    CHECK(charge2->metadata().name == std::string_view{"Charge2"});
    CHECK(charge2->metadata().full_name == std::string_view{"Charge2"});
    CHECK(charge2->metadata().publication.has_value());
    CHECK(charge2->metadata().priority == 30);

    CHECK(charge2->requirements().requires_common_parameters());
    CHECK(charge2->requirements().requires_atom_parameters());
    CHECK(!charge2->requirements().requires_bond_parameters());
    CHECK(charge2->requirements().common_parameters.size() == 6);
    CHECK(charge2->requirements().atom_parameters.size() == 3);
    CHECK(charge2->requires_parameters());
    CHECK(charge2->option_schema().size() == 1);

    CHECK(delre->id() == std::string_view{"delre"});
    CHECK(delre->metadata().name == std::string_view{"DelRe"});
    CHECK(delre->metadata().full_name == std::string_view{"Method of Del Re"});
    CHECK(delre->metadata().publication.has_value());
    CHECK(delre->metadata().priority == 130);

    CHECK(delre->requirements().requires_atom_parameters());
    CHECK(delre->requirements().requires_bond_parameters());
    CHECK(!delre->requirements().requires_common_parameters());
    CHECK(delre->requirements().atom_parameters.size() == 1);
    CHECK(delre->requirements().bond_parameters.size() == 3);
    CHECK(delre->requirements().resources.time == methods::ComplexityTerm::atoms_cubed);
    CHECK(delre->requirements().resources.memory == methods::ComplexityTerm::atoms_squared);
    CHECK(delre->requires_parameters());
    CHECK(delre->option_schema().empty());

    CHECK(mgc->id() == std::string_view{"mgc"});
    CHECK(mgc->metadata().name == std::string_view{"MGC"});
    CHECK(mgc->metadata().full_name == std::string_view{"Molecular Graph Charge"});
    CHECK(mgc->metadata().publication.has_value());
    CHECK(mgc->metadata().priority == 70);

    CHECK(!mgc->requires_parameters());
    CHECK(mgc->option_schema().empty());
    CHECK(mgc->requirements().resources.time == methods::ComplexityTerm::atoms_cubed);
    CHECK(mgc->requirements().resources.memory == methods::ComplexityTerm::atoms_squared);

    CHECK(denr->id() == std::string_view{"denr"});
    CHECK(denr->metadata().name == std::string_view{"DENR"});
    CHECK(denr->metadata().full_name == std::string_view{"Dynamical Electronegativity Relaxation"});
    CHECK(denr->metadata().publication.has_value());
    CHECK(denr->metadata().priority == 50);

    CHECK(!denr->requirements().requires_common_parameters());
    CHECK(denr->requirements().requires_atom_parameters());
    CHECK(!denr->requirements().requires_bond_parameters());
    CHECK(denr->requirements().common_parameters.empty());
    CHECK(denr->requirements().atom_parameters.size() == 2);
    CHECK(denr->requires_parameters());
    const auto denr_option_schema = denr->option_schema();
    CHECK(denr_option_schema.size() == 2);
    CHECK(denr_option_schema[0].id == std::string_view{"step"});
    CHECK(denr_option_schema[0].type == methods::MethodOptionType::floating_point);
    CHECK(std::get_if<double>(&denr_option_schema[0].default_value) != nullptr);
    CHECK(*std::get_if<double>(&denr_option_schema[0].default_value) == 0.1);
    CHECK(denr_option_schema[1].id == std::string_view{"iterations"});
    CHECK(denr_option_schema[1].type == methods::MethodOptionType::integer);
    CHECK(std::get_if<int>(&denr_option_schema[1].default_value) != nullptr);
    CHECK(*std::get_if<int>(&denr_option_schema[1].default_value) == 3);

    CHECK(denr->requirements().resources.time == methods::ComplexityTerm::atoms_cubed);
    CHECK(denr->requirements().resources.memory == methods::ComplexityTerm::atoms_squared);

    CHECK(kcm->id() == std::string_view{"kcm"});
    CHECK(kcm->metadata().name == std::string_view{"KCM"});
    CHECK(kcm->metadata().full_name == std::string_view{"Kirchhoff Charge Model"});
    CHECK(kcm->metadata().publication.has_value());
    CHECK(kcm->metadata().priority == 60);

    CHECK(kcm->requirements().requires_atom_parameters());
    CHECK(!kcm->requirements().requires_common_parameters());
    CHECK(!kcm->requirements().requires_bond_parameters());
    CHECK(kcm->requirements().atom_parameters.size() == 2);
    CHECK(kcm->requires_parameters());
    CHECK(kcm->option_schema().empty());

    CHECK(kcm->requirements().resources.time == methods::ComplexityTerm::atoms_cubed);
    CHECK(kcm->requirements().resources.memory == methods::ComplexityTerm::atoms_squared);

    CHECK(tsef->id() == std::string_view{"tsef"});
    CHECK(tsef->metadata().name == std::string_view{"TSEF"});
    CHECK(tsef->metadata().full_name ==
          std::string_view{"Topologically Symmetrical Energy Function"});
    CHECK(tsef->metadata().publication.has_value());
    CHECK(tsef->metadata().priority == 55);

    CHECK(!tsef->requirements().requires_common_parameters());
    CHECK(tsef->requirements().requires_atom_parameters());
    CHECK(!tsef->requirements().requires_bond_parameters());
    CHECK(tsef->requirements().atom_parameters.size() == 2);
    CHECK(tsef->requires_parameters());
    CHECK(tsef->option_schema().empty());

    CHECK(tsef->requirements().resources.time == methods::ComplexityTerm::atoms_cubed);
    CHECK(tsef->requirements().resources.memory == methods::ComplexityTerm::atoms_squared);

    CHECK(qeq->id() == std::string_view{"qeq"});
    CHECK(qeq->metadata().name == std::string_view{"QEq"});
    CHECK(qeq->metadata().full_name == std::string_view{"Charge Equilibration"});
    CHECK(qeq->metadata().publication.has_value());
    CHECK(qeq->metadata().priority == 170);

    CHECK(qeq->requirements().coordinates);
    CHECK(qeq->requirements().requires_atom_parameters());
    CHECK(!qeq->requirements().requires_common_parameters());
    CHECK(!qeq->requirements().requires_bond_parameters());
    CHECK(qeq->requirements().atom_parameters.size() == 2);
    CHECK(qeq->requires_parameters());
    CHECK(qeq->option_schema().size() == 1);

    CHECK(qeq->requirements().resources.time == methods::ComplexityTerm::atoms_cubed);
    CHECK(qeq->requirements().resources.memory == methods::ComplexityTerm::atoms_squared);
    CHECK(qeq->requirements().resources.supports_cutoff);
    CHECK(qeq->requirements().resources.supports_cover);
    CHECK(qeq->requirements().resources.fragment_target_charge_policy ==
          methods::FragmentTargetChargePolicy::proportional_to_atom_count);

    CHECK(sqe->id() == std::string_view{"sqe"});
    CHECK(sqe->metadata().name == std::string_view{"SQE"});
    CHECK(sqe->metadata().full_name == std::string_view{"Split-charge Equilibration"});
    CHECK(sqe->metadata().publication.has_value());
    CHECK(sqe->metadata().priority == 90);
    CHECK(sqe->requirements().coordinates);
    CHECK(sqe->requirements().requires_atom_parameters());
    CHECK(sqe->requirements().requires_bond_parameters());
    CHECK(!sqe->requirements().requires_common_parameters());
    CHECK(sqe->requirements().atom_parameters.size() == 3);
    CHECK(sqe->requirements().bond_parameters.size() == 1);
    CHECK(sqe->requires_parameters());
    CHECK(sqe->option_schema().empty());
    CHECK(sqe->requirements().resources.time == methods::ComplexityTerm::bonds_cubed);
    CHECK(sqe->requirements().resources.memory == methods::ComplexityTerm::bonds_squared);
    CHECK(sqe->requirements().resources.supports_cutoff);
    CHECK(sqe->requirements().resources.supports_cover);
    CHECK(sqe->requirements().resources.fragment_target_charge_policy ==
          methods::FragmentTargetChargePolicy::zero);

    CHECK(sqeq0->id() == std::string_view{"sqeq0"});
    CHECK(sqeq0->metadata().name == std::string_view{"SQE+q0"});
    CHECK(sqeq0->metadata().full_name ==
          std::string_view{"Split-charge Equilibration with Initial Formal Charges"});
    CHECK(sqeq0->metadata().publication.has_value());
    CHECK(sqeq0->metadata().priority == 80);
    CHECK(sqeq0->requirements().coordinates);
    CHECK(sqeq0->requirements().requires_atom_parameters());
    CHECK(sqeq0->requirements().requires_bond_parameters());
    CHECK(!sqeq0->requirements().requires_common_parameters());
    CHECK(sqeq0->requirements().atom_parameters.size() == 3);
    CHECK(sqeq0->requirements().bond_parameters.size() == 1);
    CHECK(sqeq0->requires_parameters());
    CHECK(sqeq0->option_schema().empty());
    CHECK(sqeq0->requirements().resources.time == methods::ComplexityTerm::bonds_cubed);
    CHECK(sqeq0->requirements().resources.memory == methods::ComplexityTerm::bonds_squared);
    CHECK(sqeq0->requirements().resources.supports_cutoff);
    CHECK(sqeq0->requirements().resources.supports_cover);
    CHECK(sqeq0->requirements().resources.fragment_target_charge_policy ==
          methods::FragmentTargetChargePolicy::proportional_to_atom_count);

    CHECK(sqeqp->id() == std::string_view{"sqeqp"});
    CHECK(sqeqp->metadata().name == std::string_view{"SQE+qp"});
    CHECK(sqeqp->metadata().full_name ==
          std::string_view{"Split-charge Equilibration with Parameterized Initial Charges"});
    CHECK(sqeqp->metadata().publication.has_value());
    CHECK(sqeqp->metadata().priority == 210);
    CHECK(sqeqp->requirements().coordinates);
    CHECK(sqeqp->requirements().requires_atom_parameters());
    CHECK(sqeqp->requirements().requires_bond_parameters());
    CHECK(!sqeqp->requirements().requires_common_parameters());
    CHECK(sqeqp->requirements().atom_parameters.size() == 4);
    CHECK(sqeqp->requirements().bond_parameters.size() == 1);
    CHECK(sqeqp->requires_parameters());
    CHECK(sqeqp->option_schema().empty());
    CHECK(sqeqp->requirements().resources.time == methods::ComplexityTerm::bonds_cubed);
    CHECK(sqeqp->requirements().resources.memory == methods::ComplexityTerm::bonds_squared);
    CHECK(sqeqp->requirements().resources.supports_cutoff);
    CHECK(sqeqp->requirements().resources.supports_cover);
    CHECK(sqeqp->requirements().resources.fragment_target_charge_policy ==
          methods::FragmentTargetChargePolicy::proportional_to_atom_count);

    CHECK(eem->id() == std::string_view{"eem"});
    CHECK(eem->metadata().name == std::string_view{"EEM"});
    CHECK(eem->metadata().full_name == std::string_view{"Electronegativity Equalization Method"});
    CHECK(eem->metadata().publication.has_value());
    CHECK(eem->metadata().priority == 200);

    CHECK(eem->requirements().coordinates);
    CHECK(eem->requirements().requires_common_parameters());
    CHECK(eem->requirements().requires_atom_parameters());
    CHECK(!eem->requirements().requires_bond_parameters());
    CHECK(eem->requirements().common_parameters.size() == 1);
    CHECK(eem->requirements().atom_parameters.size() == 2);
    CHECK(eem->requires_parameters());
    CHECK(eem->option_schema().empty());

    CHECK(eem->requirements().resources.time == methods::ComplexityTerm::atoms_cubed);
    CHECK(eem->requirements().resources.memory == methods::ComplexityTerm::atoms_squared);
    CHECK(eem->requirements().resources.supports_cutoff);
    CHECK(eem->requirements().resources.supports_cover);
    CHECK(eem->requirements().resources.fragment_target_charge_policy ==
          methods::FragmentTargetChargePolicy::proportional_to_atom_count);

    CHECK(smpqeq->id() == std::string_view{"smpqeq"});
    CHECK(smpqeq->metadata().name == std::string_view{"SMP/QEq"});
    CHECK(smpqeq->metadata().full_name ==
          std::string_view{"Self-Consistent Charge Equilibration Method"});
    CHECK(smpqeq->metadata().publication.has_value());
    CHECK(smpqeq->metadata().priority == 160);

    CHECK(smpqeq->requirements().coordinates);
    CHECK(smpqeq->requirements().requires_atom_parameters());
    CHECK(!smpqeq->requirements().requires_common_parameters());
    CHECK(!smpqeq->requirements().requires_bond_parameters());
    CHECK(smpqeq->requirements().atom_parameters.size() == 4);
    CHECK(smpqeq->requires_parameters());
    CHECK(smpqeq->option_schema().empty());

    CHECK(smpqeq->requirements().resources.time == methods::ComplexityTerm::atoms_cubed);
    CHECK(smpqeq->requirements().resources.memory == methods::ComplexityTerm::atoms_squared);

    CHECK(sfkeem->id() == std::string_view{"sfkeem"});
    CHECK(sfkeem->metadata().name == std::string_view{"SFKEEM"});
    CHECK(sfkeem->metadata().full_name ==
          std::string_view{"Selfconsistent Functional Kernel Equalized Electronegativity Method"});
    CHECK(sfkeem->metadata().publication.has_value());
    CHECK(sfkeem->metadata().priority == 180);

    CHECK(sfkeem->requirements().coordinates);
    CHECK(sfkeem->requirements().requires_common_parameters());
    CHECK(sfkeem->requirements().requires_atom_parameters());
    CHECK(!sfkeem->requirements().requires_bond_parameters());
    CHECK(sfkeem->requirements().common_parameters.size() == 1);
    CHECK(sfkeem->requirements().atom_parameters.size() == 2);
    CHECK(sfkeem->requires_parameters());
    CHECK(sfkeem->option_schema().empty());

    CHECK(sfkeem->requirements().resources.time == methods::ComplexityTerm::atoms_cubed);
    CHECK(sfkeem->requirements().resources.memory == methods::ComplexityTerm::atoms_squared);
    CHECK(!sfkeem->requirements().resources.supports_cutoff);

    CHECK(eqeq->id() == std::string_view{"eqeq"});
    CHECK(eqeq->metadata().name == std::string_view{"EQeq"});
    CHECK(eqeq->metadata().full_name == std::string_view{"Extended Charge Equilibration Method"});
    CHECK(eqeq->metadata().publication.has_value());
    CHECK(eqeq->metadata().priority == 150);

    CHECK(eqeq->requirements().coordinates);
    CHECK(!eqeq->requirements().requires_common_parameters());
    CHECK(!eqeq->requirements().requires_atom_parameters());
    CHECK(!eqeq->requirements().requires_bond_parameters());
    CHECK(!eqeq->requires_parameters());
    CHECK(eqeq->option_schema().empty());

    CHECK(eqeq->requirements().resources.time == methods::ComplexityTerm::atoms_cubed);
    CHECK(eqeq->requirements().resources.memory == methods::ComplexityTerm::atoms_squared);
    CHECK(eqeq->requirements().resources.supports_cutoff);
    CHECK(eqeq->requirements().resources.supports_cover);
    CHECK(eqeq->requirements().resources.fragment_target_charge_policy ==
          methods::FragmentTargetChargePolicy::proportional_to_atom_count);

    CHECK(eqeqc->id() == std::string_view{"eqeqc"});
    CHECK(eqeqc->metadata().name == std::string_view{"EQeq+C"});
    CHECK(eqeqc->metadata().full_name ==
          std::string_view{"Bond-Order-Corrected Extended Charge Equilibration Method"});
    CHECK(eqeqc->metadata().publication.has_value());
    CHECK(eqeqc->metadata().priority == 140);

    CHECK(eqeqc->requirements().coordinates);
    CHECK(eqeqc->requirements().requires_common_parameters());
    CHECK(eqeqc->requirements().requires_atom_parameters());
    CHECK(!eqeqc->requirements().requires_bond_parameters());
    CHECK(eqeqc->requirements().common_parameters.size() == 1);
    CHECK(eqeqc->requirements().atom_parameters.size() == 1);
    CHECK(eqeqc->requires_parameters());
    CHECK(eqeqc->option_schema().empty());

    CHECK(eqeqc->requirements().resources.time == methods::ComplexityTerm::atoms_cubed);
    CHECK(eqeqc->requirements().resources.memory == methods::ComplexityTerm::atoms_squared);
    CHECK(eqeqc->requirements().resources.supports_cutoff);
    CHECK(eqeqc->requirements().resources.supports_cover);
    CHECK(eqeqc->requirements().resources.fragment_target_charge_policy ==
          methods::FragmentTargetChargePolicy::proportional_to_atom_count);

    CHECK(abeem->id() == std::string_view{"abeem"});
    CHECK(abeem->metadata().name == std::string_view{"ABEEM"});
    CHECK(abeem->metadata().full_name ==
          std::string_view{"Atom-Bond Electronegativity Equalization Method"});
    CHECK(abeem->metadata().publication.has_value());
    CHECK(abeem->metadata().priority == 190);

    CHECK(abeem->requirements().coordinates);
    CHECK(abeem->requirements().requires_common_parameters());
    CHECK(abeem->requirements().requires_atom_parameters());
    CHECK(abeem->requirements().requires_bond_parameters());
    CHECK(abeem->requirements().common_parameters.size() == 1);
    CHECK(abeem->requirements().atom_parameters.size() == 3);
    CHECK(abeem->requirements().bond_parameters.size() == 4);
    CHECK(abeem->requires_parameters());
    CHECK(abeem->option_schema().empty());

    CHECK(abeem->requirements().resources.time == methods::ComplexityTerm::atoms_plus_bonds_cubed);
    CHECK(abeem->requirements().resources.memory ==
          methods::ComplexityTerm::atoms_plus_bonds_squared);
    CHECK(abeem->requirements().resources.supports_cutoff);
    CHECK(abeem->requirements().resources.supports_cover);
    CHECK(abeem->requirements().resources.fragment_target_charge_policy ==
          methods::FragmentTargetChargePolicy::proportional_to_atom_count);

    const auto water = chargefw::test::make_water_graph();

    const auto dummy_charges = calculate(*dummy, water);
    CHECK(dummy_charges.size() == water.atom_count());

    for (const auto charge : dummy_charges.values()) {
        CHECK(charge == 0.0);
    }

    const auto veem_charges = calculate(*veem, water);
    CHECK(veem_charges.size() == water.atom_count());

    CHECK(veem_charges[0] < 0.0);
    CHECK(veem_charges[1] > 0.0);
    CHECK(veem_charges[2] > 0.0);
    CHECK(std::abs(veem_charges[1] - veem_charges[2]) < 1.0e-12);
    CHECK(std::abs(veem_charges.total()) < 1.0e-12);

    const auto mgc_charges = calculate(*mgc, water);

    CHECK(mgc_charges.size() == water.atom_count());
    CHECK(mgc_charges[0] < 0.0);
    CHECK(mgc_charges[1] > 0.0);
    CHECK(mgc_charges[2] > 0.0);
    CHECK(std::abs(mgc_charges[1] - mgc_charges[2]) < 1.0e-12);
    CHECK(std::abs(mgc_charges.total()) < 1.0e-12);

    const auto charged_pair = chargefw::test::make_formally_charged_pair();

    const auto formal_charges = calculate(*formal, charged_pair);
    CHECK(formal_charges.size() == charged_pair.atom_count());

    CHECK(formal_charges[0] == 1.0);
    CHECK(formal_charges[1] == -1.0);
    CHECK(formal_charges.total() == 0.0);
}