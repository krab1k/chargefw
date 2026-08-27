#include <chargefw/calculation/assessment.h>
#include <chargefw/calculation/calculation.h>
#include <chargefw/calculation/execution_policy.h>
#include <chargefw/config.h>
#include <chargefw/core/atom.h>
#include <chargefw/core/bond.h>
#include <chargefw/core/conformer.h>
#include <chargefw/core/molecule.h>
#include <chargefw/core/molecule_collection.h>
#include <chargefw/core/position.h>
#include <chargefw/methods/method_applicability.h>
#include <chargefw/methods/method_options.h>
#include <chargefw/parameters/io/parameter_set_io.h>

#include <nanobind/nanobind.h>
#include <nanobind/ndarray.h>
#include <nanobind/stl/string.h>
#include <nanobind/stl/vector.h>

#include <cstddef>
#include <cstdint>
#include <filesystem>
#include <optional>
#include <stdexcept>
#include <string>
#include <string_view>
#include <unordered_map>
#include <utility>
#include <vector>

namespace nb = nanobind;

namespace {

using integer_array_1d =
    nb::ndarray<const std::int64_t, nb::ndim<1>, nb::c_contig, nb::device::cpu>;
using integer_array_2d =
    nb::ndarray<const std::int64_t, nb::ndim<2>, nb::c_contig, nb::device::cpu>;
using coordinate_array = nb::ndarray<const double, nb::ndim<3>, nb::c_contig, nb::device::cpu>;

auto sequence_string(const nb::sequence& values, const std::size_t index) -> std::string {
    return nb::cast<std::string>(values[index]);
}

auto make_molecule(integer_array_1d atomic_numbers, integer_array_1d formal_charges,
                   integer_array_2d bonds, coordinate_array coordinates, nb::sequence atom_names,
                   nb::sequence conformer_names, std::string name) -> chargefw::core::Molecule {
    const auto atom_count = static_cast<std::size_t>(atomic_numbers.shape(0));

    if (formal_charges.shape(0) != atomic_numbers.shape(0)) {
        throw std::invalid_argument{"formal charge count must match atomic number count"};
    }
    if (bonds.shape(1) != 3) {
        throw std::invalid_argument{"bonds must have three columns"};
    }
    if (coordinates.shape(1) != static_cast<std::int64_t>(atom_count) ||
        coordinates.shape(2) != 3) {
        throw std::invalid_argument{"coordinate shape does not match the molecule"};
    }

    std::vector<chargefw::core::Atom> atoms;
    atoms.reserve(atom_count);
    for (std::size_t index = 0; index < atom_count; ++index) {
        atoms.emplace_back(static_cast<int>(atomic_numbers.data()[index]),
                           static_cast<int>(formal_charges.data()[index]),
                           sequence_string(atom_names, index));
    }

    const auto bond_count = static_cast<std::size_t>(bonds.shape(0));
    std::vector<chargefw::core::Bond> native_bonds;
    native_bonds.reserve(bond_count);
    for (std::size_t index = 0; index < bond_count; ++index) {
        const auto offset = index * 3;
        native_bonds.emplace_back(
            static_cast<std::size_t>(bonds.data()[offset]),
            static_cast<std::size_t>(bonds.data()[offset + 1]),
            chargefw::core::bond_order_from_value(static_cast<int>(bonds.data()[offset + 2])));
    }

    const auto conformer_count = static_cast<std::size_t>(coordinates.shape(0));
    std::vector<chargefw::core::Conformer> conformers;
    conformers.reserve(conformer_count);
    for (std::size_t conformer_index = 0; conformer_index < conformer_count; ++conformer_index) {
        std::vector<chargefw::core::Position> positions;
        positions.reserve(atom_count);
        for (std::size_t atom_index = 0; atom_index < atom_count; ++atom_index) {
            const auto offset = (conformer_index * atom_count + atom_index) * 3;
            positions.push_back(chargefw::core::Position{
                .x = coordinates.data()[offset],
                .y = coordinates.data()[offset + 1],
                .z = coordinates.data()[offset + 2],
            });
        }
        conformers.emplace_back(std::move(positions),
                                sequence_string(conformer_names, conformer_index));
    }

    return chargefw::core::Molecule{std::move(atoms), std::move(native_bonds),
                                    std::move(conformers), std::move(name)};
}

auto make_collection(nb::sequence molecules, std::string name)
    -> chargefw::core::MoleculeCollection {
    std::vector<chargefw::core::Molecule> native_molecules;
    native_molecules.reserve(nb::len(molecules));
    for (std::size_t index = 0; index < nb::len(molecules); ++index) {
        native_molecules.emplace_back(nb::cast<const chargefw::core::Molecule&>(molecules[index]));
    }
    return chargefw::core::MoleculeCollection{std::move(native_molecules), std::move(name)};
}

class NativeParameterCatalog {
  public:
    explicit NativeParameterCatalog(std::vector<chargefw::parameters::ParameterSet> parameter_sets)
        : parameter_sets_{std::move(parameter_sets)} {}

    [[nodiscard]] auto size() const noexcept -> std::size_t {
        return parameter_sets_.size();
    }

    [[nodiscard]] auto parameter_sets() const noexcept
        -> const std::vector<chargefw::parameters::ParameterSet>& {
        return parameter_sets_;
    }

  private:
    std::vector<chargefw::parameters::ParameterSet> parameter_sets_;
};

auto make_parameter_catalog(const std::string& directory) -> NativeParameterCatalog {
    return NativeParameterCatalog{
        chargefw::parameters::load_parameter_sets_json_directory(std::filesystem::path{directory})};
}

auto optional_string(const nb::handle value) -> std::optional<std::string> {
    if (value.is_none()) {
        return std::nullopt;
    }
    return nb::cast<std::string>(value);
}

auto optional_double(const nb::handle value) -> std::optional<double> {
    if (value.is_none()) {
        return std::nullopt;
    }
    return nb::cast<double>(value);
}

auto optional_size(const nb::handle value) -> std::optional<std::size_t> {
    if (value.is_none()) {
        return std::nullopt;
    }
    return nb::cast<std::size_t>(value);
}

auto method_options(const nb::dict& values)
    -> std::unordered_map<std::string, chargefw::methods::MethodOptions> {
    std::unordered_map<std::string, chargefw::methods::MethodOptions> result;
    for (const auto& [method_value, options_value] : values) {
        const auto method_id = nb::cast<std::string>(method_value);
        const auto options = nb::cast<nb::dict>(options_value);
        std::unordered_map<std::string, chargefw::methods::MethodOptionValue> native_values;
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
        result.emplace(method_id, chargefw::methods::MethodOptions{std::move(native_values)});
    }
    return result;
}

auto issue_kind(const chargefw::methods::PrerequisiteIssueKind kind) -> std::string_view {
    using Kind = chargefw::methods::PrerequisiteIssueKind;
    switch (kind) {
    case Kind::invalid_options:
        return "invalid_options";
    case Kind::missing_feature:
        return "missing_feature";
    case Kind::invalid_geometry:
        return "invalid_geometry";
    case Kind::unsupported_molecule:
        return "unsupported_molecule";
    case Kind::missing_parameters:
        return "missing_parameters";
    case Kind::parameter_classification_failed:
        return "parameter_classification_failed";
    }
    return "unknown";
}

auto execution_issue_kind(const chargefw::methods::ExecutionIssueKind kind) -> std::string_view {
    using Kind = chargefw::methods::ExecutionIssueKind;
    switch (kind) {
    case Kind::resource_threshold_exceeded:
        return "resource_threshold_exceeded";
    case Kind::unsupported_execution_mode:
        return "unsupported_execution_mode";
    }
    return "unknown";
}

auto execution_availability(const chargefw::methods::ExecutionAvailability availability)
    -> std::string_view {
    using Availability = chargefw::methods::ExecutionAvailability;
    switch (availability) {
    case Availability::available:
        return "available";
    case Availability::available_with_warning:
        return "available_with_warning";
    case Availability::unsupported:
        return "unsupported";
    }
    return "unknown";
}

auto prerequisite_issue(const chargefw::methods::PrerequisiteIssue& issue) -> nb::dict {
    auto result = nb::dict{};
    result["kind"] = std::string{issue_kind(issue.kind)};
    result["message"] = issue.message;
    result["molecule_index"] =
        issue.molecule_index.has_value() ? nb::cast(*issue.molecule_index) : nb::none();
    result["atom_index"] = issue.atom_index.has_value() ? nb::cast(*issue.atom_index) : nb::none();
    result["bond_index"] = issue.bond_index.has_value() ? nb::cast(*issue.bond_index) : nb::none();
    result["conformer_index"] =
        issue.conformer_index.has_value() ? nb::cast(*issue.conformer_index) : nb::none();
    return result;
}

auto execution_issue(const chargefw::methods::ExecutionIssue& issue) -> nb::dict {
    auto result = nb::dict{};
    result["kind"] = std::string{execution_issue_kind(issue.kind)};
    result["message"] = issue.message;
    result["molecule_index"] =
        issue.molecule_index.has_value() ? nb::cast(*issue.molecule_index) : nb::none();
    return result;
}

auto execution_assessment(const chargefw::methods::ExecutionAssessment& assessment) -> nb::dict {
    auto result = nb::dict{};
    result["mode"] = std::string{chargefw::calculation::to_string(assessment.mode)};
    result["availability"] = std::string{execution_availability(assessment.availability)};
    auto issues = nb::list{};
    for (const auto& issue : assessment.issues) {
        issues.append(execution_issue(issue));
    }
    result["issues"] = std::move(issues);
    return result;
}

auto applicability_report(const chargefw::calculation::ApplicabilityReport& report) -> nb::dict {
    auto result = nb::dict{};
    auto applicable = nb::list{};
    for (const auto& candidate : report.applicable) {
        auto value = nb::dict{};
        value["method_id"] = candidate.method_id;
        value["parameter_set_id"] = candidate.parameter_set_id.has_value()
                                        ? nb::cast(std::string{*candidate.parameter_set_id})
                                        : nb::none();
        auto assessments = nb::list{};
        for (const auto& assessment : candidate.execution_assessments) {
            assessments.append(execution_assessment(assessment));
        }
        value["execution_assessments"] = std::move(assessments);
        applicable.append(std::move(value));
    }
    result["applicable"] = std::move(applicable);

    auto rejected = nb::list{};
    for (const auto& candidate : report.rejected) {
        auto value = nb::dict{};
        value["method_id"] = candidate.method_id;
        value["parameter_set_id"] = candidate.parameter_set_id.has_value()
                                        ? nb::cast(std::string{*candidate.parameter_set_id})
                                        : nb::none();
        auto issues = nb::list{};
        for (const auto& issue : candidate.issues) {
            issues.append(prerequisite_issue(issue));
        }
        value["issues"] = std::move(issues);
        rejected.append(std::move(value));
    }
    result["rejected"] = std::move(rejected);
    result["selected_candidate_index"] = report.selected_candidate_index.has_value()
                                             ? nb::cast(*report.selected_candidate_index)
                                             : nb::none();
    return result;
}

auto execution_policy(const chargefw::calculation::ExecutionPolicy& policy) -> nb::dict {
    auto result = nb::dict{};
    result["mode"] = std::string{chargefw::calculation::to_string(policy.mode())};
    result["radius"] = policy.radius().has_value() ? nb::cast(*policy.radius()) : nb::none();
    result["charge_correction"] =
        std::string{chargefw::calculation::to_string(policy.charge_correction())};
    return result;
}

auto method_option_value(const chargefw::methods::MethodOptionValue& value) -> nb::object {
    return std::visit([](const auto& item) { return nb::cast(item); }, value);
}

auto method_options(const chargefw::methods::MethodOptions& options) -> nb::dict {
    auto result = nb::dict{};
    for (const auto& [id, value] : options.values()) {
        result[id.c_str()] = method_option_value(value);
    }
    return result;
}

auto effective_calculation(const chargefw::calculation::EffectiveCalculation& effective)
    -> nb::dict {
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

auto status(const chargefw::calculation::ExecutionStatus value) -> std::string_view {
    using Status = chargefw::calculation::ExecutionStatus;
    switch (value) {
    case Status::success:
        return "success";
    case Status::invalid_input_or_request:
        return "invalid_input_or_request";
    case Status::no_executable_plan:
        return "no_executable_plan";
    case Status::numerical_failure:
        return "numerical_failure";
    case Status::cancelled:
        return "cancelled";
    }
    return "unknown";
}

auto execution_result(const chargefw::calculation::ExecutionResult& value) -> nb::dict {
    auto result = nb::dict{};
    result["status"] = std::string{status(value.status)};
    result["applicability"] = applicability_report(value.applicability);
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
    charge_set["method_id"] = std::string{charges.method_id()};
    const auto parameter_set_id = charges.parameter_set_id();
    charge_set["parameter_set_id"] =
        parameter_set_id.has_value() ? nb::cast(std::string{*parameter_set_id}) : nb::none();
    auto assignments = nb::list{};
    for (const auto& assignment : charges.assignments()) {
        auto item = nb::dict{};
        item["molecule_index"] = assignment.target.molecule_index;
        item["conformer_index"] = assignment.target.conformer_index.has_value()
                                      ? nb::cast(*assignment.target.conformer_index)
                                      : nb::none();
        item["values"] = std::vector<double>{assignment.charges.values().begin(),
                                             assignment.charges.values().end()};
        assignments.append(std::move(item));
    }
    charge_set["assignments"] = std::move(assignments);
    result["charges"] = std::move(charge_set);
    return result;
}

class NativeAssessment {
  public:
    NativeAssessment(chargefw::calculation::AssessmentResult assessment,
                     const std::size_t max_threads)
        : assessment_{std::move(assessment)}, max_threads_{max_threads} {}

    [[nodiscard]] auto report() const -> nb::dict {
        auto result = nb::dict{};
        result["applicability"] = applicability_report(assessment_.applicability());
        if (assessment_.execution_policy().has_value()) {
            result["execution_policy"] = execution_policy(*assessment_.execution_policy());
        } else {
            result["execution_policy"] = nb::none();
        }
        auto issues = nb::list{};
        for (const auto& issue : assessment_.execution_issues()) {
            issues.append(execution_issue(issue));
        }
        result["execution_issues"] = std::move(issues);
        result["applicability_seconds"] = assessment_.applicability_seconds();
        result["executable"] = assessment_.executable();
        return result;
    }

    auto calculate() -> nb::dict {
        if (consumed_) {
            throw std::runtime_error{"assessment can only be calculated once"};
        }
        consumed_ = true;
        const auto result = [&] {
            nb::gil_scoped_release release;
            return chargefw::calculation::calculate(std::move(assessment_), max_threads_);
        }();
        return execution_result(result);
    }

  private:
    chargefw::calculation::AssessmentResult assessment_;
    std::size_t max_threads_ = 0;
    bool consumed_ = false;
};

auto make_assessment(const chargefw::core::MoleculeCollection& molecules,
                     const NativeParameterCatalog& catalog, const nb::handle method_id,
                     const nb::handle parameter_set_id, const nb::dict& options,
                     const bool permissive_types, const std::string& execution,
                     const nb::handle radius, const nb::handle charge_correction,
                     const nb::handle cutoff_threshold, const nb::handle cover_threshold,
                     const std::size_t max_threads) -> NativeAssessment {
    const auto selection_kind =
        chargefw::calculation::execution_selection_kind_from_string(execution);
    const auto correction =
        charge_correction.is_none()
            ? std::optional<chargefw::calculation::ChargeCorrectionPolicy>{}
            : std::optional{chargefw::calculation::charge_correction_policy_from_string(
                  nb::cast<std::string>(charge_correction))};
    auto request = chargefw::calculation::AssessmentRequest{
        .molecules = molecules,
        .parameter_sets = catalog.parameter_sets(),
        .method_id = optional_string(method_id),
        .parameter_set_id = optional_string(parameter_set_id),
        .method_options = method_options(options),
        .classification_options = {.permissive_types = permissive_types},
        .execution_selection =
            chargefw::calculation::ExecutionSelection{selection_kind, optional_double(radius),
                                                      correction},
        .resource_policy = {.cutoff_atom_threshold = optional_size(cutoff_threshold),
                            .cover_atom_threshold = optional_size(cover_threshold),
                            .max_threads = max_threads},
    };
    nb::gil_scoped_release release;
    return NativeAssessment{chargefw::calculation::assess(std::move(request)), max_threads};
}

auto make_assessment_args(nb::args args) -> NativeAssessment {
    if (args.size() != 12) {
        throw std::invalid_argument{
            "internal assessment bridge received an unexpected argument count"};
    }
    return make_assessment(nb::cast<const chargefw::core::MoleculeCollection&>(args[0]),
                           nb::cast<const NativeParameterCatalog&>(args[1]), args[2], args[3],
                           nb::cast<nb::dict>(args[4]), nb::cast<bool>(args[5]),
                           nb::cast<std::string>(args[6]), args[7], args[8], args[9], args[10],
                           nb::cast<std::size_t>(args[11]));
}

} // namespace

NB_MODULE(_chargefw, module) {
    module.doc() = "Private native extension for ChargeFW.";
    module.def("version", [] { return CHARGEFW_VERSION_STRING; });

    nb::class_<chargefw::core::Molecule>(module, "_NativeMolecule")
        .def_prop_ro("name", &chargefw::core::Molecule::name)
        .def_prop_ro("atom_count", &chargefw::core::Molecule::atom_count)
        .def_prop_ro("bond_count", &chargefw::core::Molecule::bond_count)
        .def_prop_ro("conformer_count", &chargefw::core::Molecule::conformer_count);
    nb::class_<chargefw::core::MoleculeCollection>(module, "_NativeMoleculeCollection")
        .def_prop_ro("name", &chargefw::core::MoleculeCollection::name)
        .def_prop_ro("size", &chargefw::core::MoleculeCollection::size);
    nb::class_<NativeParameterCatalog>(module, "_NativeParameterCatalog")
        .def_prop_ro("size", &NativeParameterCatalog::size);
    nb::class_<NativeAssessment>(module, "_NativeAssessment")
        .def("report", &NativeAssessment::report)
        .def("calculate", &NativeAssessment::calculate);

    module.def("_make_molecule", &make_molecule, nb::arg("atomic_numbers"),
               nb::arg("formal_charges"), nb::arg("bonds"), nb::arg("coordinates"),
               nb::arg("atom_names"), nb::arg("conformer_names"), nb::arg("name"));
    module.def("_make_collection", &make_collection, nb::arg("molecules"), nb::arg("name"));
    module.def("_load_parameter_catalog", &make_parameter_catalog, nb::arg("directory"));
    module.def("_make_assessment", &make_assessment_args);
}
