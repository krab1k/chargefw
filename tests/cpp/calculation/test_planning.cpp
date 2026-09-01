#include "support/test_molecules.h"

#include <chargefw/calculation/calculation.h>
#include <chargefw/core/molecule_collection.h>

#include <snitch/snitch.hpp>

#include <variant>
#include <vector>

namespace calculation = chargefw::calculation;
namespace core = chargefw::core;
namespace methods = chargefw::methods;

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
