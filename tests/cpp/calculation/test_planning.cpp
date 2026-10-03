#include "support/test_molecules.h"
#include "support/test_parameters.h"

#include <chargefw/calculation/calculation.h>
#include <chargefw/core/atom.h>
#include <chargefw/core/conformer.h>
#include <chargefw/core/molecule.h>
#include <chargefw/core/molecule_collection.h>
#include <chargefw/methods/method.h>
#include <chargefw/parameters/models/common_parameters.h>
#include <chargefw/parameters/models/parameter_set.h>
#include <chargefw/parameters/models/parameter_set_metadata.h>

#include <algorithm>
#include <snitch/snitch.hpp>

#include <limits>
#include <stdexcept>
#include <string_view>
#include <utility>
#include <variant>
#include <vector>

namespace calculation = chargefw::calculation;
namespace core = chargefw::core;
namespace methods = chargefw::methods;

namespace {

auto make_isolated_ion_pair() -> core::Molecule {
    return core::Molecule{
        std::vector{core::Atom{1}, core::Atom{12, 2}},
        {},
        {core::Conformer{{core::Position{0.0, 0.0, 0.0}, core::Position{3.0, 0.0, 0.0}}}},
        "isolated-ion-pair"};
}

auto make_bonded_pair() -> core::Molecule {
    return core::Molecule{std::vector{core::Atom{1}, core::Atom{6}},
                          std::vector{core::Bond{0, 1}},
                          {},
                          "bonded-pair"};
}

auto make_eem_parameters() -> chargefw::parameters::ParameterSet {
    return chargefw::parameters::ParameterSet{
        chargefw::parameters::ParameterSetMetadata{
            .id = "test-eem", .method_id = "eem", .name = "Test EEM"},
        chargefw::parameters::CommonParameters{{{.name = "kappa", .value = 1.0}}},
        chargefw::parameters::AtomParameters{
            {{.key = chargefw::test::plain_atom_key(1),
              .parameters = {{.name = "A", .value = 1.0}, {.name = "B", .value = 10.0}}},
             {.key = chargefw::test::plain_atom_key(8),
              .parameters = {{.name = "A", .value = 2.0}, {.name = "B", .value = 10.0}}}}}};
}

} // namespace

TEST_CASE("explicit unsupported execution has no selected plan or fallback",
          "[calculation][planning]") {
    auto assessment = calculation::assess(calculation::AssessmentRequest{
        .molecules = core::MoleculeCollection{std::vector{chargefw::test::make_water()}},
        .method_id = "formal",
        .execution_selection = calculation::ExecutionSelection{
            calculation::ExecutionSelectionKind::cover, calculation::minimum_reduced_radius}});

    CHECK(assessment.plans().empty());
    REQUIRE(assessment.rejections().size() == 1);
    REQUIRE(assessment.rejections()[0].policy.has_value());
    CHECK(assessment.rejections()[0].policy->mode() == calculation::ExecutionMode::cover);
    REQUIRE(assessment.rejections()[0].issues.size() == 1);
    CHECK(std::get<methods::ExecutionIssue>(assessment.rejections()[0].issues[0]).kind ==
          methods::ExecutionIssueKind::unsupported_execution_mode);

    const auto result = calculation::calculate(assessment);
    CHECK(result.status == calculation::ExecutionStatus::no_executable_plan);
    CHECK_FALSE(result.calculated());
}

TEST_CASE("explicit no-plan assessment reports rejected scientific prerequisites",
          "[calculation][planning]") {
    auto assessment = calculation::assess(calculation::AssessmentRequest{
        .molecules = core::MoleculeCollection{std::vector{chargefw::test::make_water()}},
        .method_id = "smpqeq"});

    CHECK(assessment.plans().empty());
    REQUIRE(assessment.rejections().size() == 1);
    CHECK_FALSE(assessment.rejections()[0].policy.has_value());
    CHECK(assessment.rejections()[0].method_id == "smpqeq");
    CHECK_FALSE(assessment.rejections()[0].issues.empty());

    const auto result = calculation::calculate(assessment);
    CHECK(result.status == calculation::ExecutionStatus::no_executable_plan);
    CHECK_FALSE(result.calculated());
}

TEST_CASE("fixed charge embedding sources fail closed before planning", "[calculation][planning]") {
    auto make_request = [](const bool explicit_method) {
        auto request = calculation::AssessmentRequest{
            .molecules = core::MoleculeCollection{std::vector{make_isolated_ion_pair()}},
            .fixed_charge_embedding = calculation::FixedChargeEmbedding{
                .sources = {{.molecule_index = 0, .atom_index = 1, .charge = 2.0}},
                .charge_provenance = "test charges"}};
        if (explicit_method) {
            request.method_id = "eem";
            request.execution_selection =
                calculation::ExecutionSelection{calculation::ExecutionSelectionKind::full};
        }
        return request;
    };
    const auto check_rejected = [](auto&& request) {
        auto rejected = false;
        try {
            static_cast<void>(calculation::assess(std::forward<decltype(request)>(request)));
        } catch (const std::invalid_argument& error) {
            rejected = true;
            CHECK(std::string_view{error.what()} ==
                  "fixed charge embedding sources are not supported by assessment planning");
        }
        CHECK(rejected);
    };

    auto explicit_full = make_request(true);
    check_rejected(explicit_full);
    REQUIRE(explicit_full.fixed_charge_embedding.has_value());
    REQUIRE(explicit_full.fixed_charge_embedding->sources.size() == 1);
    CHECK(explicit_full.fixed_charge_embedding->sources[0].charge == 2.0);
    check_rejected(make_request(false));
    check_rejected(make_request(true));
}

TEST_CASE("fixed charge embedding source selectors are validated", "[calculation][planning]") {
    const auto check_invalid = [](std::vector<core::Molecule> molecules,
                                  std::vector<calculation::FixedAtomCharge> sources,
                                  const std::string_view diagnostic) {
        auto request = calculation::AssessmentRequest{
            .molecules = core::MoleculeCollection{std::move(molecules)},
            .fixed_charge_embedding =
                calculation::FixedChargeEmbedding{.sources = std::move(sources)}};
        auto rejected = false;
        try {
            static_cast<void>(calculation::assess(std::move(request)));
        } catch (const std::invalid_argument& error) {
            rejected = true;
            CHECK(std::string_view{error.what()}.contains(diagnostic));
        }
        CHECK(rejected);
    };

    const auto ion_pair = make_isolated_ion_pair();
    check_invalid({ion_pair}, {{2, 0, 1.0}}, "molecule index 2 is out of range");
    check_invalid({ion_pair}, {{std::numeric_limits<std::size_t>::max(), 0, 1.0}},
                  "is out of range (molecule count 1)");
    check_invalid({ion_pair}, {{0, std::numeric_limits<std::size_t>::max(), 1.0}},
                  "out of range (atom count 2)");
    check_invalid({ion_pair}, {{0, 0, 1.0}, {0, 0, -1.0}}, "duplicate fixed charge source");
    for (const auto charge :
         {std::numeric_limits<double>::quiet_NaN(), std::numeric_limits<double>::infinity(),
          -std::numeric_limits<double>::infinity()}) {
        check_invalid({ion_pair}, {{0, 0, charge}}, "has non-finite charge");
    }
    check_invalid({make_bonded_pair()}, {{0, 0, 0.5}}, "is bonded to atom 1");
    check_invalid({make_bonded_pair()}, {{0, 1, 0.5}}, "is bonded to atom 0");
    check_invalid({make_bonded_pair()}, {{0, 0, 0.5}, {0, 1, -0.5}}, "is bonded to atom 1");
    check_invalid({core::Molecule{std::vector{core::Atom{6}}}, ion_pair}, {{0, 0, 0.0}},
                  "leaves no active atoms in molecule 0");
}

TEST_CASE("valid fixed charge selectors reach the unsupported planning gate",
          "[calculation][planning]") {
    auto request = calculation::AssessmentRequest{
        .molecules = core::MoleculeCollection{std::vector{
            make_isolated_ion_pair(), make_isolated_ion_pair(), make_isolated_ion_pair()}},
        .fixed_charge_embedding = calculation::FixedChargeEmbedding{
            .sources = {{.molecule_index = 1, .atom_index = 1, .charge = -0.5},
                        {.molecule_index = 0, .atom_index = 1, .charge = 0.0},
                        {.molecule_index = 2, .atom_index = 0, .charge = 0.25}}}};
    try {
        static_cast<void>(calculation::assess(std::move(request)));
        CHECK(false);
    } catch (const std::invalid_argument& error) {
        CHECK(std::string_view{error.what()} ==
              "fixed charge embedding sources are not supported by assessment planning");
    }
}

TEST_CASE("empty fixed charge embedding is equivalent to absence", "[calculation][planning]") {
    auto make_request = [](const bool with_empty_embedding) {
        auto request = calculation::AssessmentRequest{
            .molecules = core::MoleculeCollection{std::vector{chargefw::test::make_water()}},
            .parameter_sets = {make_eem_parameters()}};
        if (with_empty_embedding) {
            request.fixed_charge_embedding =
                calculation::FixedChargeEmbedding{.charge_provenance = "ignored when empty"};
        }
        return request;
    };

    const auto ordinary = calculation::assess(make_request(false));
    const auto empty = calculation::assess(make_request(true));
    REQUIRE_FALSE(ordinary.plans().empty());
    REQUIRE(empty.plans().size() == ordinary.plans().size());
    for (std::size_t index = 0; index < ordinary.plans().size(); ++index) {
        CHECK(empty.plans()[index].candidate().method->id() ==
              ordinary.plans()[index].candidate().method->id());
        CHECK(empty.plans()[index].policy().mode() == ordinary.plans()[index].policy().mode());
        CHECK(empty.plans()[index].policy().radius() == ordinary.plans()[index].policy().radius());
    }
    const auto ordinary_result = calculation::calculate(ordinary);
    const auto empty_result = calculation::calculate(empty);
    REQUIRE(ordinary_result.calculated());
    REQUIRE(empty_result.calculated());
    CHECK(std::ranges::equal(empty_result.charges->assignment(0).charges.values(),
                             ordinary_result.charges->assignment(0).charges.values()));
}
