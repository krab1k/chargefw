#include "support/test_molecules.h"

#include <chargefw/features/prepared_molecule.h>
#include <chargefw/methods/calculation_input.h>
#include <chargefw/methods/method_options.h>
#include <chargefw/methods/method_registry.h>

#include <cassert>
#include <cmath>
#include <string>
#include <string_view>
#include <vector>

namespace features = chargefw::features;
namespace methods = chargefw::methods;

namespace {

auto calculate(const methods::Method& method, const chargefw::core::Molecule& molecule)
    -> chargefw::charges::AtomicCharges {
    const features::PreparedMolecule prepared_molecule{molecule};
    const auto method_options = methods::make_default_options(method.option_schema());

    const methods::CalculationInput input{prepared_molecule, method_options};

    return method.calculate(input);
}

} // namespace

auto main() -> int {
    const auto& registry = methods::method_registry();

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
    const auto* eem = registry.find("eem");
    const auto* smpqeq = registry.find("smpqeq");

    assert(dummy != nullptr);
    assert(formal != nullptr);
    assert(veem != nullptr);
    assert(peoe != nullptr);
    assert(mpeoe != nullptr);
    assert(gdac != nullptr);
    assert(charge2 != nullptr);
    assert(delre != nullptr);
    assert(mgc != nullptr);
    assert(denr != nullptr);
    assert(kcm != nullptr);
    assert(tsef != nullptr);
    assert(qeq != nullptr);
    assert(eem != nullptr);
    assert(smpqeq != nullptr);

    const auto method_names = registry.names();

    assert((method_names == std::vector<std::string>{"charge2", "delre", "denr", "dummy", "eem",
                                                     "formal", "gdac", "kcm", "mgc", "mpeoe",
                                                     "peoe", "qeq", "smpqeq", "tsef", "veem"}));

    assert(dummy->id() == std::string_view{"dummy"});
    assert(dummy->metadata().name == std::string_view{"Dummy method"});
    assert(dummy->metadata().full_name == std::string_view{"Dummy zero charges"});
    assert(!dummy->metadata().publication.has_value());
    assert(dummy->metadata().priority == 10);
    assert(!dummy->requirements().formal_charges);
    assert(!dummy->requires_parameters());
    assert(dummy->option_schema().empty());

    assert(formal->id() == std::string_view{"formal"});
    assert(formal->metadata().name == std::string_view{"Formal"});
    assert(formal->metadata().full_name == std::string_view{"Formal atomic charges"});
    assert(!formal->metadata().publication.has_value());
    assert(formal->metadata().priority == 10);
    assert(formal->requirements().formal_charges);
    assert(!formal->requires_parameters());
    assert(formal->option_schema().empty());

    assert(veem->id() == std::string_view{"veem"});
    assert(veem->metadata().name == std::string_view{"VEEM"});
    assert(veem->metadata().full_name == std::string_view{"Valence Electrons Equalization Method"});
    assert(veem->metadata().publication.has_value());
    assert(veem->metadata().priority == 20);
    assert(veem->requirements().element_properties);
    assert(veem->requirements().resources.time == methods::ComplexityTerm::atoms);
    assert(veem->requirements().resources.memory == methods::ComplexityTerm::constant);
    assert(!veem->requires_parameters());
    assert(veem->option_schema().empty());

    assert(peoe->id() == std::string_view{"peoe"});
    assert(peoe->metadata().name == std::string_view{"PEOE"});
    assert(peoe->metadata().full_name ==
           std::string_view{"Partial Equalization of Atomic Electronegativity"});
    assert(peoe->metadata().publication.has_value());
    assert(peoe->metadata().priority == 120);

    assert(mpeoe->id() == std::string_view{"mpeoe"});
    assert(mpeoe->metadata().name == std::string_view{"MPEOE"});
    assert(mpeoe->metadata().full_name ==
           std::string_view{"Modified Partial Equalization of Atomic Electronegativity"});
    assert(mpeoe->metadata().publication.has_value());
    assert(mpeoe->metadata().priority == 110);

    assert(gdac->id() == std::string_view{"gdac"});
    assert(gdac->metadata().name == std::string_view{"GDAC"});
    assert(gdac->metadata().full_name == std::string_view{"Geometry-Dependent Net Atomic Charges"});
    assert(gdac->metadata().publication.has_value());
    assert(gdac->metadata().priority == 100);

    assert(gdac->requirements().bond_graph);
    assert(gdac->requirements().coordinates);
    assert(gdac->requirements().element_properties);
    assert(gdac->requirements().requires_atom_parameters());
    assert(!gdac->requirements().requires_common_parameters());
    assert(!gdac->requirements().requires_bond_parameters());
    assert(gdac->requirements().atom_parameters.size() == 2);
    assert(gdac->requires_parameters());
    assert(gdac->option_schema().size() == 1);

    assert(mpeoe->requirements().bond_graph);
    assert(mpeoe->requirements().requires_common_parameters());
    assert(mpeoe->requirements().requires_atom_parameters());
    assert(mpeoe->requirements().requires_bond_parameters());
    assert(mpeoe->requirements().common_parameters.size() == 1);
    assert(mpeoe->requirements().atom_parameters.size() == 2);
    assert(mpeoe->requirements().bond_parameters.size() == 1);
    assert(mpeoe->requires_parameters());
    assert(mpeoe->option_schema().size() == 1);

    assert(peoe->requirements().bond_graph);
    assert(peoe->requirements().requires_common_parameters());
    assert(peoe->requirements().requires_atom_parameters());
    assert(peoe->requirements().common_parameters.size() == 1);
    assert(peoe->requirements().atom_parameters.size() == 3);
    assert(peoe->requires_parameters());
    assert(peoe->option_schema().size() == 1);

    assert(charge2->id() == std::string_view{"charge2"});
    assert(charge2->metadata().name == std::string_view{"Charge2"});
    assert(charge2->metadata().full_name == std::string_view{"Charge2"});
    assert(charge2->metadata().publication.has_value());
    assert(charge2->metadata().priority == 30);

    assert(charge2->requirements().bond_graph);
    assert(charge2->requirements().topological_distances);
    assert(charge2->requirements().element_properties);
    assert(charge2->requirements().requires_common_parameters());
    assert(charge2->requirements().requires_atom_parameters());
    assert(!charge2->requirements().requires_bond_parameters());
    assert(charge2->requirements().common_parameters.size() == 6);
    assert(charge2->requirements().atom_parameters.size() == 3);
    assert(charge2->requires_parameters());
    assert(charge2->option_schema().size() == 1);

    assert(delre->id() == std::string_view{"delre"});
    assert(delre->metadata().name == std::string_view{"DelRe"});
    assert(delre->metadata().full_name == std::string_view{"Method of Del Re"});
    assert(delre->metadata().publication.has_value());
    assert(delre->metadata().priority == 130);

    assert(delre->requirements().bond_graph);
    assert(delre->requirements().requires_atom_parameters());
    assert(delre->requirements().requires_bond_parameters());
    assert(!delre->requirements().requires_common_parameters());
    assert(delre->requirements().atom_parameters.size() == 1);
    assert(delre->requirements().bond_parameters.size() == 3);
    assert(delre->requirements().resources.time == methods::ComplexityTerm::atoms_cubed);
    assert(delre->requirements().resources.memory == methods::ComplexityTerm::atoms_squared);
    assert(delre->requirements().resources.reject_large_without_reduction);
    assert(delre->requires_parameters());
    assert(delre->option_schema().empty());

    assert(mgc->id() == std::string_view{"mgc"});
    assert(mgc->metadata().name == std::string_view{"MGC"});
    assert(mgc->metadata().full_name == std::string_view{"Molecular Graph Charge"});
    assert(mgc->metadata().publication.has_value());
    assert(mgc->metadata().priority == 70);

    assert(mgc->requirements().bond_graph);
    assert(mgc->requirements().bond_orders);
    assert(mgc->requirements().element_properties);
    assert(!mgc->requires_parameters());
    assert(mgc->option_schema().empty());
    assert(mgc->requirements().resources.time == methods::ComplexityTerm::atoms_cubed);
    assert(mgc->requirements().resources.memory == methods::ComplexityTerm::atoms_squared);
    assert(mgc->requirements().resources.reject_large_without_reduction);

    assert(denr->id() == std::string_view{"denr"});
    assert(denr->metadata().name == std::string_view{"DENR"});
    assert(denr->metadata().full_name ==
           std::string_view{"Dynamical Electronegativity Relaxation"});
    assert(denr->metadata().publication.has_value());
    assert(denr->metadata().priority == 50);

    assert(denr->requirements().bond_graph);
    assert(denr->requirements().requires_common_parameters());
    assert(denr->requirements().requires_atom_parameters());
    assert(!denr->requirements().requires_bond_parameters());
    assert(denr->requirements().common_parameters.size() == 2);
    assert(denr->requirements().atom_parameters.size() == 2);
    assert(denr->requires_parameters());
    assert(denr->option_schema().empty());

    assert(denr->requirements().resources.time == methods::ComplexityTerm::atoms_cubed);
    assert(denr->requirements().resources.memory == methods::ComplexityTerm::atoms_squared);
    assert(denr->requirements().resources.reject_large_without_reduction);

    assert(kcm->id() == std::string_view{"kcm"});
    assert(kcm->metadata().name == std::string_view{"KCM"});
    assert(kcm->metadata().full_name == std::string_view{"Kirchhoff Charge Model"});
    assert(kcm->metadata().publication.has_value());
    assert(kcm->metadata().priority == 60);

    assert(kcm->requirements().bond_graph);
    assert(kcm->requirements().requires_atom_parameters());
    assert(!kcm->requirements().requires_common_parameters());
    assert(!kcm->requirements().requires_bond_parameters());
    assert(kcm->requirements().atom_parameters.size() == 2);
    assert(kcm->requires_parameters());
    assert(kcm->option_schema().empty());

    assert(kcm->requirements().resources.time == methods::ComplexityTerm::atoms_cubed);
    assert(kcm->requirements().resources.memory == methods::ComplexityTerm::atoms_squared);
    assert(kcm->requirements().resources.reject_large_without_reduction);

    assert(tsef->id() == std::string_view{"tsef"});
    assert(tsef->metadata().name == std::string_view{"TSEF"});
    assert(tsef->metadata().full_name ==
           std::string_view{"Topologically Symmetrical Energy Function"});
    assert(tsef->metadata().publication.has_value());
    assert(tsef->metadata().priority == 55);

    assert(tsef->requirements().bond_graph);
    assert(tsef->requirements().topological_distances);
    assert(tsef->requirements().formal_charges);
    assert(!tsef->requirements().requires_common_parameters());
    assert(tsef->requirements().requires_atom_parameters());
    assert(!tsef->requirements().requires_bond_parameters());
    assert(tsef->requirements().atom_parameters.size() == 2);
    assert(tsef->requires_parameters());
    assert(tsef->option_schema().empty());

    assert(tsef->requirements().resources.time == methods::ComplexityTerm::atoms_cubed);
    assert(tsef->requirements().resources.memory == methods::ComplexityTerm::atoms_squared);
    assert(tsef->requirements().resources.reject_large_without_reduction);

    assert(qeq->id() == std::string_view{"qeq"});
    assert(qeq->metadata().name == std::string_view{"QEq"});
    assert(qeq->metadata().full_name == std::string_view{"Charge Equilibration"});
    assert(qeq->metadata().publication.has_value());
    assert(qeq->metadata().priority == 170);

    assert(qeq->requirements().coordinates);
    assert(qeq->requirements().formal_charges);
    assert(qeq->requirements().requires_atom_parameters());
    assert(!qeq->requirements().requires_common_parameters());
    assert(!qeq->requirements().requires_bond_parameters());
    assert(qeq->requirements().atom_parameters.size() == 2);
    assert(qeq->requires_parameters());
    assert(qeq->option_schema().size() == 1);

    assert(qeq->requirements().resources.time == methods::ComplexityTerm::atoms_cubed);
    assert(qeq->requirements().resources.memory == methods::ComplexityTerm::atoms_squared);
    assert(qeq->requirements().resources.reject_large_without_reduction);

    assert(eem->id() == std::string_view{"eem"});
    assert(eem->metadata().name == std::string_view{"EEM"});
    assert(eem->metadata().full_name == std::string_view{"Electronegativity Equalization Method"});
    assert(eem->metadata().publication.has_value());
    assert(eem->metadata().priority == 200);

    assert(eem->requirements().coordinates);
    assert(eem->requirements().formal_charges);
    assert(eem->requirements().requires_common_parameters());
    assert(eem->requirements().requires_atom_parameters());
    assert(!eem->requirements().requires_bond_parameters());
    assert(eem->requirements().common_parameters.size() == 1);
    assert(eem->requirements().atom_parameters.size() == 2);
    assert(eem->requires_parameters());
    assert(eem->option_schema().empty());

    assert(eem->requirements().resources.time == methods::ComplexityTerm::atoms_cubed);
    assert(eem->requirements().resources.memory == methods::ComplexityTerm::atoms_squared);
    assert(eem->requirements().resources.reject_large_without_reduction);

    assert(smpqeq->id() == std::string_view{"smpqeq"});
    assert(smpqeq->metadata().name == std::string_view{"SMP/QEq"});
    assert(smpqeq->metadata().full_name ==
           std::string_view{"Self-Consistent Charge Equilibration Method"});
    assert(smpqeq->metadata().publication.has_value());
    assert(smpqeq->metadata().priority == 160);

    assert(smpqeq->requirements().coordinates);
    assert(smpqeq->requirements().formal_charges);
    assert(smpqeq->requirements().requires_atom_parameters());
    assert(!smpqeq->requirements().requires_common_parameters());
    assert(!smpqeq->requirements().requires_bond_parameters());
    assert(smpqeq->requirements().atom_parameters.size() == 4);
    assert(smpqeq->requires_parameters());
    assert(smpqeq->option_schema().empty());

    assert(smpqeq->requirements().resources.time == methods::ComplexityTerm::atoms_cubed);
    assert(smpqeq->requirements().resources.memory == methods::ComplexityTerm::atoms_squared);
    assert(smpqeq->requirements().resources.reject_large_without_reduction);

    const auto water = chargefw::test::make_water();

    const auto dummy_charges = calculate(*dummy, water);
    assert(dummy_charges.size() == water.atom_count());

    for (const auto charge : dummy_charges.values()) {
        assert(charge == 0.0);
    }

    const auto veem_charges = calculate(*veem, water);
    assert(veem_charges.size() == water.atom_count());

    assert(veem_charges[0] < 0.0);
    assert(veem_charges[1] > 0.0);
    assert(veem_charges[2] > 0.0);
    assert(std::abs(veem_charges[1] - veem_charges[2]) < 1.0e-12);
    assert(std::abs(veem_charges.total()) < 1.0e-12);

    const auto mgc_charges = calculate(*mgc, water);

    assert(mgc_charges.size() == water.atom_count());
    assert(mgc_charges[0] < 0.0);
    assert(mgc_charges[1] > 0.0);
    assert(mgc_charges[2] > 0.0);
    assert(std::abs(mgc_charges[1] - mgc_charges[2]) < 1.0e-12);
    assert(std::abs(mgc_charges.total()) < 1.0e-12);

    const auto charged_pair = chargefw::test::make_formally_charged_pair();

    const auto formal_charges = calculate(*formal, charged_pair);
    assert(formal_charges.size() == charged_pair.atom_count());

    assert(formal_charges[0] == 1.0);
    assert(formal_charges[1] == -1.0);
    assert(formal_charges.total() == 0.0);

    return 0;
}