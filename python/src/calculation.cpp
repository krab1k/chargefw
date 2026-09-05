#include "bindings.h"
#include "native_execution_result.h"
#include "native_parameter_catalog.h"

#include <chargefw/calculation/assessment.h>
#include <chargefw/calculation/calculation.h>
#include <chargefw/calculation/execution_policy.h>
#include <chargefw/calculation/observer.h>
#include <chargefw/core/molecule_collection.h>
#include <chargefw/methods/method.h>
#include <chargefw/methods/method_applicability.h>
#include <chargefw/methods/method_options.h>

#include <nanobind/nanobind.h>
#include <nanobind/ndarray.h>
#include <nanobind/stl/optional.h>
#include <nanobind/stl/string.h>
#include <nanobind/stl/vector.h>

#include <cstddef>
#include <memory>
#include <optional>
#include <span>
#include <stdexcept>
#include <string>
#include <string_view>
#include <unordered_map>
#include <utility>
#include <variant>
#include <vector>

namespace nb = nanobind;

namespace chargefw::python {
namespace {

auto method_options(const nb::dict& values)
    -> std::unordered_map<std::string, methods::MethodOptions> {
    std::unordered_map<std::string, methods::MethodOptions> result;
    for (const auto& [method_value, options_value] : values) {
        const auto method_id = nb::cast<std::string>(method_value);
        const auto options = nb::cast<nb::dict>(options_value);
        std::unordered_map<std::string, methods::MethodOptionValue> native_values;
        for (const auto& [option_value, value] : options) {
            const auto option_id = nb::cast<std::string>(option_value);
            if (nb::isinstance<nb::bool_>(value)) {
                native_values.emplace(option_id, nb::cast<bool>(value));
            } else if (nb::isinstance<nb::int_>(value)) {
                native_values.emplace(option_id, nb::cast<int>(value));
            } else if (nb::isinstance<nb::float_>(value)) {
                native_values.emplace(option_id, nb::cast<double>(value));
            } else if (nb::isinstance<nb::str>(value)) {
                native_values.emplace(option_id, nb::cast<std::string>(value));
            } else {
                throw std::invalid_argument{
                    "method option values must be bool, int, float, or str"};
            }
        }
        result.emplace(method_id, methods::MethodOptions{std::move(native_values)});
    }
    return result;
}

auto charge_correction_selection(const std::optional<std::string>& value)
    -> std::optional<calculation::ChargeCorrectionPolicy> {
    if (!value.has_value()) {
        return std::nullopt;
    }
    return calculation::charge_correction_policy_from_string(*value);
}

auto prerequisite_issue(const methods::PrerequisiteIssue& issue) -> nb::dict {
    auto result = nb::dict{};
    result["kind"] = std::string{methods::to_string(issue.kind)};
    result["message"] = issue.message;
    result["molecule_index"] =
        issue.molecule_index.has_value() ? nb::cast(*issue.molecule_index) : nb::none();
    result["atom_index"] = issue.atom_index.has_value() ? nb::cast(*issue.atom_index) : nb::none();
    result["bond_index"] = issue.bond_index.has_value() ? nb::cast(*issue.bond_index) : nb::none();
    result["conformer_index"] =
        issue.conformer_index.has_value() ? nb::cast(*issue.conformer_index) : nb::none();
    return result;
}

auto execution_issue(const methods::ExecutionIssue& issue) -> nb::dict {
    auto result = nb::dict{};
    result["kind"] = std::string{methods::to_string(issue.kind)};
    result["message"] = issue.message;
    result["molecule_index"] =
        issue.molecule_index.has_value() ? nb::cast(*issue.molecule_index) : nb::none();
    return result;
}

auto execution_policy(const calculation::ExecutionPolicy& policy) -> nb::dict {
    auto result = nb::dict{};
    result["mode"] = std::string{calculation::to_string(policy.mode())};
    result["radius"] = policy.radius().has_value() ? nb::cast(*policy.radius()) : nb::none();
    result["charge_correction"] = std::string{calculation::to_string(policy.charge_correction())};
    return result;
}

auto method_option_value(const methods::MethodOptionValue& value) -> nb::object {
    return std::visit([](const auto& item) { return nb::cast(item); }, value);
}

auto method_options(const methods::MethodOptions& options) -> nb::dict {
    auto result = nb::dict{};
    for (const auto& [id, value] : options.values()) {
        result[id.c_str()] = method_option_value(value);
    }
    return result;
}

auto effective_calculation(const calculation::EffectiveCalculation& effective) -> nb::dict {
    auto result = nb::dict{};
    result["method_id"] = effective.method_id;
    result["parameter_set_id"] = effective.parameter_set_id.has_value()
                                     ? nb::cast(std::string{*effective.parameter_set_id})
                                     : nb::none();
    result["method_options"] = method_options(effective.method_options);
    result["execution_policy"] = execution_policy(effective.execution_policy);
    auto issues = nb::list{};
    for (const auto& issue : effective.execution_issues) {
        issues.append(execution_issue(issue));
    }
    result["execution_issues"] = std::move(issues);
    return result;
}

auto numpy_charge_values(const std::span<const double> values)
    -> nb::ndarray<nb::numpy, const double, nb::shape<-1>> {
    auto owned_values = std::make_unique<std::vector<double>>(values.begin(), values.end());
    auto owner = nb::capsule{owned_values.get(), [](void* pointer) noexcept {
                                 delete static_cast<std::vector<double>*>(pointer);
                             }};
    auto* released_values = owned_values.release();
    return {released_values->data(), {released_values->size()}, owner};
}

auto rejection(const calculation::Rejection& rejected) -> nb::dict {
    auto result = nb::dict{};
    result["method_id"] = rejected.method_id;
    result["parameter_set_id"] =
        rejected.parameter_set_id.has_value() ? nb::cast(*rejected.parameter_set_id) : nb::none();
    result["execution_policy"] =
        rejected.policy.has_value() ? nb::cast(execution_policy(*rejected.policy)) : nb::none();
    auto prerequisite_issues = nb::list{};
    auto execution_issues = nb::list{};
    for (const auto& issue : rejected.issues) {
        if (const auto* prerequisite = std::get_if<methods::PrerequisiteIssue>(&issue)) {
            prerequisite_issues.append(prerequisite_issue(*prerequisite));
        } else {
            execution_issues.append(execution_issue(std::get<methods::ExecutionIssue>(issue)));
        }
    }
    result["prerequisite_issues"] = std::move(prerequisite_issues);
    result["execution_issues"] = std::move(execution_issues);
    return result;
}

auto execution_result(const calculation::ExecutionResult& value) -> nb::dict {
    auto result = nb::dict{};
    result["status"] = std::string{calculation::to_string(value.status)};
    auto rejections = nb::list{};
    for (const auto& rejected : value.rejections) {
        rejections.append(rejection(rejected));
    }
    result["rejections"] = std::move(rejections);
    result["failure_message"] =
        value.failure_message.has_value() ? nb::cast(*value.failure_message) : nb::none();
    auto metrics = nb::dict{};
    metrics["applicability_seconds"] = value.metrics.applicability_seconds;
    metrics["computation_seconds"] = value.metrics.computation_seconds;
    result["metrics"] = std::move(metrics);
    if (value.effective.has_value()) {
        result["effective"] = effective_calculation(*value.effective);
    } else {
        result["effective"] = nb::none();
    }

    if (!value.charges.has_value()) {
        result["charges"] = nb::none();
        return result;
    }

    const auto& charges = *value.charges;
    auto charge_set = nb::dict{};
    auto assignments = nb::list{};
    for (const auto& assignment : charges.assignments()) {
        auto item = nb::dict{};
        item["molecule_index"] = assignment.target.molecule_index;
        item["conformer_index"] = assignment.target.conformer_index.has_value()
                                      ? nb::cast(*assignment.target.conformer_index)
                                      : nb::none();
        item["values"] = numpy_charge_values(assignment.charges.values());
        assignments.append(std::move(item));
    }
    charge_set["assignments"] = std::move(assignments);
    result["charges"] = std::move(charge_set);
    return result;
}

auto calculation_phase(const calculation::CalculationPhase phase) -> std::string_view {
    switch (phase) {
    case calculation::CalculationPhase::computation_started:
        return "computation_started";
    case calculation::CalculationPhase::computation_finished:
        return "computation_finished";
    case calculation::CalculationPhase::target_started:
        return "target_started";
    case calculation::CalculationPhase::target_finished:
        return "target_finished";
    case calculation::CalculationPhase::fragment_progress:
        return "fragment_progress";
    }
    throw std::invalid_argument{"unknown calculation phase"};
}

auto calculation_progress(const calculation::CalculationProgress& progress) -> nb::dict {
    auto result = nb::dict{};
    result["phase"] = std::string{calculation_phase(progress.phase)};
    result["mode"] = std::string{calculation::to_string(progress.mode)};
    result["method_id"] = std::string{progress.method_id};
    result["target_index"] = progress.target_index;
    result["target_count"] = progress.target_count;
    result["completed_fragment_count"] = progress.completed_fragment_count;
    result["fragment_count"] = progress.fragment_count;
    result["molecule_index"] = progress.molecule_index;
    result["conformer_index"] =
        progress.conformer_index.has_value() ? nb::cast(*progress.conformer_index) : nb::none();
    result["elapsed_seconds"] = progress.elapsed_seconds;
    return result;
}

class PythonCalculationObserver final : public calculation::CalculationObserver {
  public:
    explicit PythonCalculationObserver(nb::object observer) : observer_{std::move(observer)} {}

    void on_progress(const calculation::CalculationProgress& progress) const override {
        const auto acquire = nb::gil_scoped_acquire{};
        try {
            observer_.attr("on_progress")(calculation_progress(progress));
        } catch (nb::python_error& error) {
            error.discard_as_unraisable(observer_);
        }
    }

    [[nodiscard]] auto cancelled() const noexcept -> bool override {
        const auto acquire = nb::gil_scoped_acquire{};
        try {
            return nb::cast<bool>(observer_.attr("cancelled")());
        } catch (nb::python_error& error) {
            error.discard_as_unraisable(observer_);
        } catch (...) { // NOLINT(bugprone-empty-catch): cancellation callbacks cannot throw.
        }
        return false;
    }

  private:
    nb::object observer_;
};

class NativeAssessmentState {
  public:
    NativeAssessmentState(calculation::AssessmentResult assessment, const std::size_t max_threads)
        : assessment_{std::move(assessment)}, max_threads_{max_threads} {}

    calculation::AssessmentResult assessment_;
    std::size_t max_threads_ = 0;
};

auto execution_plan(const calculation::ExecutionPlan& plan) -> nb::dict {
    const auto& candidate = plan.candidate();
    auto result = nb::dict{};
    result["method_id"] = std::string{candidate.method->id()};
    result["parameter_set_id"] = candidate.parameter_set == nullptr
                                     ? nb::none()
                                     : nb::cast(std::string{candidate.parameter_set->id()});
    result["method_options"] = method_options(candidate.method_options);
    result["execution_policy"] = execution_policy(plan.policy());
    auto warnings = nb::list{};
    for (const auto& warning : plan.warnings()) {
        warnings.append(execution_issue(warning));
    }
    result["warnings"] = std::move(warnings);
    return result;
}

class NativePlan {
  public:
    NativePlan(std::shared_ptr<NativeAssessmentState> state, const std::size_t index)
        : state_{std::move(state)}, index_{index} {}

    [[nodiscard]] auto report() const -> nb::dict {
        return execution_plan(state_->assessment_.plans()[index_]);
    }

    [[nodiscard]] auto calculate(const std::optional<std::size_t> max_threads,
                                 nb::object observer) const -> NativeExecutionResult {
        auto python_observer = std::optional<PythonCalculationObserver>{};
        const auto* native_observer = std::addressof(calculation::default_calculation_observer());
        if (!observer.is_none()) {
            python_observer.emplace(std::move(observer));
            native_observer = std::addressof(*python_observer);
        }
        auto result = [&] {
            nb::gil_scoped_release release;
            return calculation::calculate(state_->assessment_, state_->assessment_.plans()[index_],
                                          max_threads.value_or(state_->max_threads_),
                                          *native_observer);
        }();
        return NativeExecutionResult{std::move(result)};
    }

  private:
    std::shared_ptr<NativeAssessmentState> state_;
    std::size_t index_ = 0;
};

class NativeAssessment {
  public:
    NativeAssessment(calculation::AssessmentResult assessment, const std::size_t max_threads)
        : state_{std::make_shared<NativeAssessmentState>(std::move(assessment), max_threads)} {}

    [[nodiscard]] auto report() const -> nb::dict {
        auto result = nb::dict{};
        auto rejections = nb::list{};
        for (const auto& rejected : state_->assessment_.rejections()) {
            rejections.append(rejection(rejected));
        }
        result["rejections"] = std::move(rejections);
        result["applicability_seconds"] = state_->assessment_.applicability_seconds();
        return result;
    }

    [[nodiscard]] auto plans() const -> nb::list {
        auto result = nb::list{};
        for (std::size_t index = 0; index < state_->assessment_.plans().size(); ++index) {
            result.append(NativePlan{state_, index});
        }
        return result;
    }

    [[nodiscard]] auto calculate_default(nb::object observer) const -> NativeExecutionResult {
        auto python_observer = std::optional<PythonCalculationObserver>{};
        const auto* native_observer = std::addressof(calculation::default_calculation_observer());
        if (!observer.is_none()) {
            python_observer.emplace(std::move(observer));
            native_observer = std::addressof(*python_observer);
        }
        auto result = [&] {
            nb::gil_scoped_release release;
            return calculation::calculate(state_->assessment_, state_->max_threads_,
                                          *native_observer);
        }();
        return NativeExecutionResult{std::move(result)};
    }

  private:
    std::shared_ptr<NativeAssessmentState> state_;
};

auto make_assessment(const nb::sequence& molecules, std::string molecule_collection_name,
                     const NativeParameterCatalog& catalog, std::optional<std::string> method_id,
                     std::optional<std::string> parameter_set_id, const nb::dict& options,
                     const bool permissive_types, const std::string& execution,
                     const std::optional<double> radius,
                     const std::optional<std::string>& charge_correction,
                     const std::optional<std::size_t> cutoff_threshold,
                     const std::optional<std::size_t> cover_threshold,
                     const std::size_t max_threads) -> NativeAssessment {
    // Convert Python values to stable native references before releasing the GIL. The remaining
    // work copies native-owned input values and prepares the assessment without accessing Python
    // objects.
    const auto molecule_count = static_cast<std::size_t>(nb::len(molecules));
    auto source_molecules = std::vector<const core::Molecule*>{};
    source_molecules.reserve(molecule_count);
    for (std::size_t index = 0; index < molecule_count; ++index) {
        source_molecules.push_back(&nb::cast<const core::Molecule&>(molecules[index]));
    }
    auto native_method_options = method_options(options);
    nb::gil_scoped_release release;
    auto owned_molecules = std::vector<core::Molecule>{};
    owned_molecules.reserve(source_molecules.size());
    for (const auto* molecule : source_molecules) {
        owned_molecules.push_back(*molecule);
    }
    auto request = calculation::AssessmentRequest{
        .molecules = core::MoleculeCollection{std::move(owned_molecules),
                                              std::move(molecule_collection_name)},
        .parameter_sets = catalog.parameter_sets(),
        .method_id = std::move(method_id),
        .parameter_set_id = std::move(parameter_set_id),
        .method_options = std::move(native_method_options),
        .classification_options = {.permissive_types = permissive_types},
        .execution_selection =
            calculation::ExecutionSelection{
                calculation::execution_selection_kind_from_string(execution), radius,
                charge_correction_selection(charge_correction)},
        .resource_policy = {.cutoff_atom_threshold = cutoff_threshold,
                            .cover_atom_threshold = cover_threshold,
                            .max_threads = max_threads},
    };
    return NativeAssessment{calculation::assess(std::move(request)), max_threads};
}

} // namespace

void bind_calculation(nb::module_& module) {
    nb::class_<NativeExecutionResult>(module, "_NativeExecutionResult")
        .def("report",
             [](const NativeExecutionResult& value) { return execution_result(value.result()); });
    nb::class_<NativePlan>(module, "_NativePlan")
        .def("report", &NativePlan::report)
        .def("calculate", &NativePlan::calculate, nb::arg("max_threads") = nb::none(),
             nb::arg("observer") = nb::none());
    nb::class_<NativeAssessment>(module, "_NativeAssessment")
        .def("report", &NativeAssessment::report)
        .def("plans", &NativeAssessment::plans)
        .def("calculate_default", &NativeAssessment::calculate_default,
             nb::arg("observer") = nb::none());

    module.def("_make_assessment", &make_assessment, nb::arg("molecules"),
               nb::arg("molecule_collection_name"), nb::arg("catalog"), nb::arg("method_id"),
               nb::arg("parameter_set_id"), nb::arg("method_options"), nb::arg("permissive_types"),
               nb::arg("execution"), nb::arg("radius"), nb::arg("charge_correction"),
               nb::arg("cutoff_threshold"), nb::arg("cover_threshold"), nb::arg("max_threads"));
}

} // namespace chargefw::python
