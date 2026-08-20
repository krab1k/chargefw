#include "support/test_molecules.h"
#include "support/test_parameters.h"

#include <chargefw/charges/atomic_charges.h>
#include <chargefw/core/atom.h>
#include <chargefw/core/bond.h>
#include <chargefw/core/molecule.h>
#include <chargefw/core/molecule_collection.h>
#include <chargefw/features/prepared_molecule_collection.h>
#include <chargefw/methods/method.h>
#include <chargefw/methods/method_applicability.h>
#include <chargefw/methods/method_options.h>
#include <chargefw/methods/method_registry.h>
#include <chargefw/parameters/models/atom_parameters.h>
#include <chargefw/parameters/models/parameter_key.h>
#include <chargefw/parameters/models/parameter_set.h>
#include <chargefw/parameters/models/parameter_set_metadata.h>

#include <cassert>
#include <optional>
#include <span>
#include <stdexcept>
#include <string>
#include <string_view>
#include <utility>
#include <vector>

namespace core = chargefw::core;
namespace calculation = chargefw::calculation;
namespace features = chargefw::features;
namespace methods = chargefw::methods;
namespace parameters = chargefw::parameters;

namespace {

class AtomParameterMethod final : public methods::Method {
  public:
    [[nodiscard]] auto metadata() const noexcept -> const methods::MethodMetadata& override {
        static constexpr methods::MethodMetadata metadata{.id = "atom-parameter-test",
                                                          .name = "Atom parameter test",
                                                          .full_name = "Atom parameter test",
                                                          .publication = std::nullopt,
                                                          .priority = 0};

        return metadata;
    }

    [[nodiscard]] auto requirements() const -> methods::MethodRequirements override {
        auto requirements = methods::MethodRequirements{};
        requirements.bond_graph = true;
        requirements.atom_parameters = {"value"};
        return requirements;
    }

    [[nodiscard]] auto option_schema() const noexcept
        -> std::span<const methods::MethodOptionSpec> override {
        return {};
    }

    [[nodiscard]] auto calculate(const methods::CalculationInput& /* unused */) const
        -> chargefw::charges::AtomicCharges override {
        return chargefw::charges::AtomicCharges{std::vector<double>{}};
    }
};

class ResourceMethod final : public methods::Method {
  public:
    explicit ResourceMethod(methods::ResourceRequirements resources = {},
                            const bool requires_coordinates = false)
        : resources_{resources}, requires_coordinates_{requires_coordinates} {}

    [[nodiscard]] auto metadata() const noexcept -> const methods::MethodMetadata& override {
        static constexpr methods::MethodMetadata metadata{.id = "resource-test",
                                                          .name = "Resource test",
                                                          .full_name = "Resource test",
                                                          .publication = std::nullopt,
                                                          .priority = 0};
        return metadata;
    }

    [[nodiscard]] auto requirements() const -> methods::MethodRequirements override {
        auto requirements = methods::MethodRequirements{};
        requirements.coordinates = requires_coordinates_;
        requirements.resources = resources_;
        return requirements;
    }

    [[nodiscard]] auto option_schema() const noexcept
        -> std::span<const methods::MethodOptionSpec> override {
        return {};
    }

    [[nodiscard]] auto calculate(const methods::CalculationInput& /* unused */) const
        -> chargefw::charges::AtomicCharges override {
        return chargefw::charges::AtomicCharges{std::vector<double>{}};
    }

  private:
    methods::ResourceRequirements resources_;
    bool requires_coordinates_ = false;
};

auto assessment_for(const methods::ApplicableMethod& candidate,
                    const calculation::ExecutionMode mode) -> const methods::ExecutionAssessment& {
    for (const auto& assessment : candidate.execution_assessments) {
        if (assessment.mode == mode) {
            return assessment;
        }
    }

    throw std::runtime_error{"missing execution assessment"};
}

auto make_collection() -> core::MoleculeCollection {
    std::vector molecules{chargefw::test::make_water(),
                          chargefw::test::make_formally_charged_pair()};

    return core::MoleculeCollection{std::move(molecules), "test-collection"};
}

auto make_wrong_method_parameters() -> parameters::ParameterSet {
    return parameters::ParameterSet{
        parameters::ParameterSetMetadata{.id = "wrong-method-parameters",
                                         .method_id = "other-method",
                                         .name = "Wrong method parameters"},
        {},
        parameters::AtomParameters{
            {{.key = chargefw::test::atom_key(1, parameters::AtomParameterClassificationKind::PLAIN,
                                              "*"),
              .parameters = {{.name = "value", .value = 1.0}}}}}};
}

auto make_water_only_parameters() -> parameters::ParameterSet {
    return parameters::ParameterSet{
        parameters::ParameterSetMetadata{.id = "water-only-parameters",
                                         .method_id = "atom-parameter-test",
                                         .name = "Water only parameters"},
        {},
        parameters::AtomParameters{
            {{.key = chargefw::test::atom_key(
                  1, parameters::AtomParameterClassificationKind::BONDED_ELEMENTS, "O"),
              .parameters = {{.name = "value", .value = 1.0}}},
             {.key = chargefw::test::atom_key(
                  8, parameters::AtomParameterClassificationKind::BONDED_ELEMENTS, "HH"),
              .parameters = {{.name = "value", .value = 2.0}}}}}};
}

auto make_collection_parameters() -> parameters::ParameterSet {
    return parameters::ParameterSet{
        parameters::ParameterSetMetadata{.id = "collection-parameters",
                                         .method_id = "atom-parameter-test",
                                         .name = "Collection parameters"},
        {},
        parameters::AtomParameters{
            {{.key = chargefw::test::atom_key(
                  1, parameters::AtomParameterClassificationKind::BONDED_ELEMENTS, "O"),
              .parameters = {{.name = "value", .value = 1.0}}},
             {.key = chargefw::test::atom_key(
                  8, parameters::AtomParameterClassificationKind::BONDED_ELEMENTS, "HH"),
              .parameters = {{.name = "value", .value = 2.0}}},
             {.key = chargefw::test::atom_key(7, parameters::AtomParameterClassificationKind::PLAIN,
                                              "*"),
              .parameters = {{.name = "value", .value = 3.0}}},
             {.key = chargefw::test::atom_key(
                  17, parameters::AtomParameterClassificationKind::PLAIN, "*"),
              .parameters = {{.name = "value", .value = 4.0}}}}}};
}

auto make_double_bonded_carbons() -> core::Molecule {
    return core::Molecule{std::vector{core::Atom{6}, core::Atom{6}},
                          std::vector{core::Bond{0, 1, core::BondOrder::DOUBLE}},
                          {},
                          "double-bonded-carbons"};
}

auto make_permissive_parameters(const bool include_exact_match) -> parameters::ParameterSet {
    auto entries = std::vector<parameters::AtomParameterEntry>{};

    if (include_exact_match) {
        entries.push_back(
            {.key = chargefw::test::atom_key(
                 6, parameters::AtomParameterClassificationKind::HIGHEST_BOND_ORDER, "2"),
             .parameters = {{.name = "value", .value = 2.0}}});
    }

    entries.push_back({.key = chargefw::test::atom_key(
                           6, parameters::AtomParameterClassificationKind::HIGHEST_BOND_ORDER, "1"),
                       .parameters = {{.name = "value", .value = 1.0}}});

    return parameters::ParameterSet{
        parameters::ParameterSetMetadata{.id = "permissive-parameters",
                                         .method_id = "atom-parameter-test",
                                         .name = "Permissive parameters"},
        {},
        parameters::AtomParameters{std::move(entries)}};
}

} // namespace

auto main() -> int {
    const auto collection = make_collection();
    const features::PreparedMoleculeCollection prepared_collection{collection};

    const auto& registry = methods::method_registry();
    const auto* dummy = registry.find("dummy");
    const auto* mgc = registry.find("mgc");

    assert(dummy != nullptr);
    assert(mgc != nullptr);

    const AtomParameterMethod atom_parameter_method;

    const std::vector<const methods::Method*> candidate_methods{dummy, &atom_parameter_method};

    const std::vector parameter_sets{make_wrong_method_parameters(), make_water_only_parameters(),
                                     make_collection_parameters()};

    const auto result = methods::find_applicable_methods({.molecules = prepared_collection,
                                                          .methods = candidate_methods,
                                                          .parameter_sets = parameter_sets});

    assert(!result.empty());
    assert(result.applicable.size() == 2);

    const auto& dummy_applicable = result.applicable[0];

    assert(dummy_applicable.method == dummy);
    assert(dummy_applicable.parameter_set == nullptr);
    assert(!dummy_applicable.uses_parameters());
    assert(dummy_applicable.classifications.empty());
    assert(dummy_applicable.execution_assessments.size() == 3);
    assert(assessment_for(dummy_applicable, calculation::ExecutionMode::full).availability ==
           methods::ExecutionAvailability::available);
    assert(assessment_for(dummy_applicable, calculation::ExecutionMode::cutoff).availability ==
           methods::ExecutionAvailability::unsupported);
    assert(assessment_for(dummy_applicable, calculation::ExecutionMode::cover).availability ==
           methods::ExecutionAvailability::unsupported);

    const auto& parameterized_applicable = result.applicable[1];

    assert(parameterized_applicable.method == &atom_parameter_method);
    assert(parameterized_applicable.parameter_set != nullptr);
    assert(parameterized_applicable.uses_parameters());
    assert(parameterized_applicable.parameter_set->id() ==
           std::string_view{"collection-parameters"});
    assert(parameterized_applicable.classifications.size() == prepared_collection.size());

    const auto& water_classification = parameterized_applicable.classifications[0];

    assert(water_classification.atom().size() == 3);
    assert(water_classification.atom()[0] == 1);
    assert(water_classification.atom()[1] == 0);
    assert(water_classification.atom()[2] == 0);

    const auto& charged_pair_classification = parameterized_applicable.classifications[1];

    assert(charged_pair_classification.atom().size() == 2);
    assert(charged_pair_classification.atom()[0] == 2);
    assert(charged_pair_classification.atom()[1] == 3);

    assert(result.rejected.size() == 1);

    assert(result.rejected[0].method_index == 1);
    assert(result.rejected[0].parameter_set_index == std::optional<std::size_t>{1});
    assert(!result.rejected[0].issues.empty());
    assert(result.rejected[0].issues[0].kind ==
           methods::PrerequisiteIssueKind::parameter_classification_failed);

    const core::MoleculeCollection double_bonded_collection{
        std::vector{make_double_bonded_carbons()}};
    const features::PreparedMoleculeCollection double_bonded_prepared{double_bonded_collection};
    const std::vector<const methods::Method*> permissive_methods{&atom_parameter_method};
    const std::vector permissive_only_parameters{make_permissive_parameters(false)};

    const auto strict_result =
        methods::find_applicable_methods({.molecules = double_bonded_prepared,
                                          .methods = permissive_methods,
                                          .parameter_sets = permissive_only_parameters});
    assert(strict_result.empty());
    assert(strict_result.rejected.size() == 1);
    assert(strict_result.rejected[0].issues[0].kind ==
           methods::PrerequisiteIssueKind::parameter_classification_failed);

    const auto permissive_result =
        methods::find_applicable_methods({.molecules = double_bonded_prepared,
                                          .methods = permissive_methods,
                                          .parameter_sets = permissive_only_parameters,
                                          .classification_options = {.permissive_types = true}});
    assert(permissive_result.applicable.size() == 1);
    assert(permissive_result.applicable[0].classifications[0].atom()[0] == 0);
    assert(permissive_result.applicable[0].classifications[0].atom()[1] == 0);
    assert(assessment_for(permissive_result.applicable[0], calculation::ExecutionMode::full)
               .availability == methods::ExecutionAvailability::available);

    const std::vector exact_and_permissive_parameters{make_permissive_parameters(true)};
    const auto exact_result =
        methods::find_applicable_methods({.molecules = double_bonded_prepared,
                                          .methods = permissive_methods,
                                          .parameter_sets = exact_and_permissive_parameters,
                                          .classification_options = {.permissive_types = true}});
    assert(exact_result.applicable.size() == 1);
    assert(exact_result.applicable[0].classifications[0].atom()[0] == 0);
    assert(exact_result.applicable[0].classifications[0].atom()[1] == 0);
    assert(
        assessment_for(exact_result.applicable[0], calculation::ExecutionMode::full).availability ==
        assessment_for(permissive_result.applicable[0], calculation::ExecutionMode::full)
            .availability);

    const ResourceMethod inexpensive_method{};
    const std::vector<const methods::Method*> inexpensive_methods{&inexpensive_method};
    const auto inexpensive_result =
        methods::find_applicable_methods({.molecules = prepared_collection,
                                          .methods = inexpensive_methods,
                                          .parameter_sets = {},
                                          .resource_policy = {.full_atom_threshold = 1}});
    assert(inexpensive_result.applicable.size() == 1);
    assert(assessment_for(inexpensive_result.applicable[0], calculation::ExecutionMode::full)
               .availability == methods::ExecutionAvailability::available);

    const ResourceMethod expensive_method{
        methods::ResourceRequirements{.time = methods::ComplexityTerm::atoms_cubed,
                                      .memory = methods::ComplexityTerm::atoms_squared}};
    const std::vector<const methods::Method*> expensive_methods{&expensive_method};
    const auto below_threshold_result =
        methods::find_applicable_methods({.molecules = prepared_collection,
                                          .methods = expensive_methods,
                                          .parameter_sets = {},
                                          .resource_policy = {.full_atom_threshold = 3}});
    const auto& below_threshold_full =
        assessment_for(below_threshold_result.applicable[0], calculation::ExecutionMode::full);
    assert(below_threshold_full.availability == methods::ExecutionAvailability::available);
    assert(below_threshold_full.issues.empty());

    const auto threshold_result =
        methods::find_applicable_methods({.molecules = prepared_collection,
                                          .methods = expensive_methods,
                                          .parameter_sets = {},
                                          .resource_policy = {.full_atom_threshold = 2}});
    assert(threshold_result.applicable.size() == 1);
    const auto& threshold_full =
        assessment_for(threshold_result.applicable[0], calculation::ExecutionMode::full);
    assert(threshold_full.availability == methods::ExecutionAvailability::available_with_warning);
    assert(threshold_full.issues.size() == 1);
    assert(threshold_full.issues[0].kind ==
           methods::ExecutionIssueKind::resource_threshold_exceeded);
    assert(threshold_full.issues[0].molecule_index == std::optional<std::size_t>{0});
    assert(assessment_for(threshold_result.applicable[0], calculation::ExecutionMode::cutoff)
               .availability == methods::ExecutionAvailability::unsupported);
    assert(assessment_for(threshold_result.applicable[0], calculation::ExecutionMode::cover)
               .availability == methods::ExecutionAvailability::unsupported);

    const std::vector<const methods::Method*> topology_methods{mgc};
    const auto topology_result =
        methods::find_applicable_methods({.molecules = prepared_collection,
                                          .methods = topology_methods,
                                          .parameter_sets = {},
                                          .resource_policy = {.full_atom_threshold = 2}});
    assert(topology_result.applicable.size() == 1);
    assert(assessment_for(topology_result.applicable[0], calculation::ExecutionMode::full)
               .availability == methods::ExecutionAvailability::available_with_warning);
    assert(assessment_for(topology_result.applicable[0], calculation::ExecutionMode::cutoff)
               .availability == methods::ExecutionAvailability::unsupported);
    assert(assessment_for(topology_result.applicable[0], calculation::ExecutionMode::cover)
               .availability == methods::ExecutionAvailability::unsupported);

    const auto collection_threshold_result =
        methods::find_applicable_methods({.molecules = prepared_collection,
                                          .methods = expensive_methods,
                                          .parameter_sets = {},
                                          .resource_policy = {.full_atom_threshold = 1}});
    const auto& collection_threshold_full =
        assessment_for(collection_threshold_result.applicable[0], calculation::ExecutionMode::full);
    assert(collection_threshold_full.issues.size() == prepared_collection.size());

    const auto unlimited_result = methods::find_applicable_methods(
        {.molecules = prepared_collection,
         .methods = expensive_methods,
         .parameter_sets = {},
         .resource_policy = {.full_atom_threshold = std::nullopt}});
    const auto& unlimited_full =
        assessment_for(unlimited_result.applicable[0], calculation::ExecutionMode::full);
    assert(unlimited_full.availability == methods::ExecutionAvailability::available);
    assert(unlimited_full.issues.empty());

    const ResourceMethod topology_reduced_method{
        methods::ResourceRequirements{.supports_cutoff = true, .supports_cover = true}};
    const std::vector<const methods::Method*> topology_reduced_methods{&topology_reduced_method};
    const auto topology_reduced_result =
        methods::find_applicable_methods({.molecules = prepared_collection,
                                          .methods = topology_reduced_methods,
                                          .parameter_sets = {}});
    assert(assessment_for(topology_reduced_result.applicable[0], calculation::ExecutionMode::cutoff)
               .availability == methods::ExecutionAvailability::unsupported);
    assert(assessment_for(topology_reduced_result.applicable[0], calculation::ExecutionMode::cover)
               .availability == methods::ExecutionAvailability::unsupported);

    const ResourceMethod spatial_reduced_method{
        {.supports_cutoff = true,
         .supports_cover = true,
         .fragment_target_charge_policy =
             methods::FragmentTargetChargePolicy::proportional_to_atom_count},
        true};
    const std::vector<const methods::Method*> spatial_reduced_methods{&spatial_reduced_method};
    const core::MoleculeCollection spatial_collection{std::vector{chargefw::test::make_water()}};
    const features::PreparedMoleculeCollection spatial_prepared_collection{spatial_collection};
    const auto spatial_reduced_result =
        methods::find_applicable_methods({.molecules = spatial_prepared_collection,
                                          .methods = spatial_reduced_methods,
                                          .parameter_sets = {}});
    assert(assessment_for(spatial_reduced_result.applicable[0], calculation::ExecutionMode::cutoff)
               .availability == methods::ExecutionAvailability::available);
    assert(assessment_for(spatial_reduced_result.applicable[0], calculation::ExecutionMode::cover)
               .availability == methods::ExecutionAvailability::available);

    return 0;
}
