#include "support/test_molecules.h"

#include <chargefw/adapters/generated_output.h>
#include <chargefw/adapters/molecule_record.h>
#include <chargefw/charges/atomic_charges.h>
#include <chargefw/charges/charge_collection.h>
#include <snitch/snitch.hpp>

#include <sstream>
#include <string>
#include <vector>

namespace generated_output = chargefw::adapters::generated_output;
namespace adapters = chargefw::adapters;
namespace charges = chargefw::charges;

namespace {

auto records() -> std::vector<adapters::ImportedMoleculeRecord> {
    return {{.molecule = chargefw::test::make_two_conformer_water(),
             .identity = {.record_id = "water"}}};
}

auto invariant_charges() -> charges::ChargeSet {
    return charges::ChargeSet{
        "formal",
        {{.target = {.molecule_index = 0}, .charges = charges::AtomicCharges{{-0.8, 0.4, 0.4}}}}};
}

} // namespace

TEST_CASE("generated molecular output applies shared selection policy", "[adapters][output]") {
    const auto molecules = records();
    const auto charge_set = invariant_charges();

    std::ostringstream sdf;
    generated_output::write(sdf, molecules, charge_set, generated_output::Format::sdf_v3000,
                            "ChargeFW", "test");
    CHECK(sdf.str().contains("V3000"));
    CHECK(sdf.str().contains("M  V30 2 H 0.757 0.5859 0 0"));
    CHECK_FALSE(sdf.str().contains("M  V30 2 H 1.1 0 0 0"));
    CHECK(sdf.str().contains("CHARGEFW_CHARGES_1"));

    std::ostringstream mol2;
    generated_output::write(mol2, molecules, charge_set, generated_output::Format::mol2);
    CHECK(mol2.str().contains("USER_CHARGES"));

    std::ostringstream mmcif;
    generated_output::write(mmcif, molecules, charge_set, generated_output::Format::mmcif);
    CHECK(mmcif.str().contains("_sb_ncbr_partial_atomic_charges."));
}
