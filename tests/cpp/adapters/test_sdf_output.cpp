#include <chargefw/adapters/native/sdf_output.h>
#include <chargefw/charges/atomic_charges.h>
#include <chargefw/charges/charge_collection.h>
#include <chargefw/core/atom.h>
#include <chargefw/core/bond.h>
#include <chargefw/core/conformer.h>
#include <chargefw/core/molecule.h>
#include <chargefw/core/position.h>
#include <snitch/snitch.hpp>

#include <array>
#include <sstream>
#include <stdexcept>
#include <string>
#include <vector>

namespace sdf_output = chargefw::adapters::native::sdf_output;
namespace charges = chargefw::charges;

namespace {

[[nodiscard]] auto assignment(const std::size_t molecule_index, std::vector<double> values)
    -> charges::ChargeAssignment {
    return {.target = {.molecule_index = molecule_index, .conformer_index = 0},
            .charges = charges::AtomicCharges{std::move(values)}};
}

constexpr std::string_view first_record =
    "first\nchargefw\n\n  2  0  0  0  0  0  0  0  0  0999 V2000\n"
    "    0.0000    0.0000    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n"
    "    1.0000    0.0000    0.0000 H   0  0  0  0  0  0  0  0  0  0  0  0\n"
    "M  END\n> <SOURCE_ID>\nkept\n\n> <CHARGEFW_CHARGES_1>\n9.0000 9.0000\n\n";

} // namespace

TEST_CASE("SDF output preserves source records while replacing charge properties",
          "[adapters][sdf]") {
    const auto second_record =
        std::string{"second\r\nchargefw\r\n\r\n  1  0  0  0  0  0  0  0  0  0999 V2000\r\n"
                    "    0.0000    0.0000    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0\r\n"
                    "M  END\r\n"};
    const auto source = std::string{first_record} + "$$$$\n" + second_record + "$$$$\r\n";
    const auto type_one = std::vector{assignment(0, {-0.12345, 0.12345}), assignment(1, {-0.5})};
    const auto type_two = std::vector{assignment(0, {-0.2, 0.2}), assignment(1, {-0.4})};
    const auto properties =
        std::array{sdf_output::ChargeProperty{.charge_type_id = 1,
                                              .assignments = type_one,
                                              .method = "eem",
                                              .parameter_set = "2015",
                                              .software_name = "ChargeFW",
                                              .software_version = "0.0.1"},
                   sdf_output::ChargeProperty{.charge_type_id = 2, .assignments = type_two}};

    {
        auto output = std::ostringstream{};
        sdf_output::SdfWriter{output}.write_preserving_buffer(source, properties,
                                                              sdf_output::WriteMode::replace);
        const auto text = output.str();
        CHECK(text.contains("> <SOURCE_ID>\nkept\n\n"));
        CHECK_FALSE(text.contains("9.0000 9.0000"));
        CHECK(text.contains("> <CHARGEFW_CHARGES_1>\n-0.1235 0.1235\n\n"));
        CHECK(text.contains("> <CHARGEFW_CHARGES_2>\n-0.2000 0.2000\n\n"));
        CHECK(text.contains("> <CHARGEFW_CHARGE_METADATA_1>\ntype=empirical; method=eem; "
                            "parameter_set=2015; software_name=ChargeFW; "
                            "software_version=0.0.1\n\n"));
        CHECK(text.contains("> <CHARGEFW_CHARGES_1>\r\n-0.5000\r\n\r\n"
                            "> <CHARGEFW_CHARGE_METADATA_1>\r\n"
                            "type=empirical; method=eem; parameter_set=2015; "
                            "software_name=ChargeFW; software_version=0.0.1\r\n\r\n"
                            "> <CHARGEFW_CHARGES_2>\r\n-0.4000\r\n\r\n"
                            "> <CHARGEFW_CHARGE_METADATA_2>\r\n"
                            "type=empirical; method=unknown; parameter_set=.; "
                            "software_name=unknown; software_version=unknown\r\n\r\n$$$$\r\n"));
        auto expected_prefix = std::string{first_record};
        const auto old_field = std::string{"> <CHARGEFW_CHARGES_1>\n9.0000 9.0000\n\n"};
        expected_prefix.erase(expected_prefix.find(old_field), old_field.size());
        CHECK(text.starts_with(expected_prefix));
    }

    {
        const auto source_with_metadata = std::string{first_record} +
                                          "> <CHARGEFW_CHARGE_METADATA_1>\ntype=empirical; "
                                          "method=obsolete\n\n$$$$\n";
        auto output = std::ostringstream{};
        const auto one_record_assignments = std::vector{assignment(0, {-0.12345, 0.12345})};
        const auto one_record =
            std::array{sdf_output::ChargeProperty{.charge_type_id = 1,
                                                  .assignments = one_record_assignments,
                                                  .method = "eem",
                                                  .parameter_set = "2015",
                                                  .software_name = "ChargeFW",
                                                  .software_version = "0.0.1"}};
        sdf_output::SdfWriter{output}.write_preserving_buffer(source_with_metadata, one_record,
                                                              sdf_output::WriteMode::replace);
        const auto text = output.str();
        CHECK_FALSE(text.contains("method=obsolete"));
        CHECK(text.contains("method=eem; parameter_set=2015; software_name=ChargeFW; "
                            "software_version=0.0.1"));
    }

    {
        auto output = std::ostringstream{};
        const auto only_type_one =
            std::array{sdf_output::ChargeProperty{.charge_type_id = 1, .assignments = type_one}};
        sdf_output::SdfWriter{output}.write_preserving_buffer(source, only_type_one,
                                                              sdf_output::WriteMode::append);
        const auto text = output.str();
        CHECK(text.contains("> <CHARGEFW_CHARGES_1>\n9.0000 9.0000\n\n"));
        CHECK(text.contains("> <CHARGEFW_CHARGES_1>\n-0.1235 0.1235\n\n"));
        CHECK(
            text.contains("> <CHARGEFW_CHARGE_METADATA_1>\ntype=empirical; method=unknown; "
                          "parameter_set=.; software_name=unknown; software_version=unknown\n\n"));
    }

    {
        auto rejected = false;
        try {
            auto output = std::ostringstream{};
            const auto invalid = std::array{
                sdf_output::ChargeProperty{.charge_type_id = 0, .assignments = type_one}};
            sdf_output::SdfWriter{output}.write_preserving_buffer(source, invalid);
        } catch (const std::invalid_argument&) {
            rejected = true;
        }
        CHECK(rejected);
    }

    const auto molecule = chargefw::core::Molecule{
        {chargefw::core::Atom{8, -1}, chargefw::core::Atom{1}, chargefw::core::Atom{1}},
        {chargefw::core::Bond{0, 1}, chargefw::core::Bond{0, 2}},
        {chargefw::core::Conformer{{{.x = 0.0, .y = 0.0, .z = 0.0},
                                    {.x = 0.9, .y = 0.0, .z = 0.0},
                                    {.x = -0.2, .y = 0.9, .z = 0.0}}}},
        "hydroxide"};
    const auto generated_type_one = std::vector{assignment(0, {-0.8, 0.4, 0.4})};
    const auto generated_type_two = std::vector{assignment(0, {-0.7, 0.35, 0.35})};
    const auto generated_properties = std::array{
        sdf_output::ChargeProperty{.charge_type_id = 1, .assignments = generated_type_one},
        sdf_output::ChargeProperty{.charge_type_id = 2, .assignments = generated_type_two}};

    {
        auto output = std::ostringstream{};
        sdf_output::SdfWriter{output}.write_generated(molecule, generated_properties,
                                                      sdf_output::MolFormat::v2000);
        const auto text = output.str();
        CHECK(text.contains("  3  2  0  0  0  0  0  0  0  0999 V2000\n"));
        CHECK(text.contains("M  CHG  1   1  -1\nM  END\n"));
        CHECK(text.contains("> <CHARGEFW_CHARGES_1>\n-0.8000 0.4000 0.4000\n\n"));
        CHECK(
            text.contains("> <CHARGEFW_CHARGE_METADATA_1>\ntype=empirical; method=unknown; "
                          "parameter_set=.; software_name=unknown; software_version=unknown\n\n"));
        CHECK(text.ends_with("$$$$\n"));
    }

    {
        auto output = std::ostringstream{};
        sdf_output::SdfWriter{output}.write_generated(molecule, generated_properties,
                                                      sdf_output::MolFormat::v3000);
        const auto text = output.str();
        CHECK(text.contains("999 V3000\nM  V30 BEGIN CTAB\nM  V30 COUNTS 3 2 0 0 0\n"));
        CHECK(text.contains("M  V30 1 O 0 0 0 0 CHG=-1\n"));
        CHECK(text.contains("M  V30 1 1 1 2\nM  V30 2 1 1 3\n"));
        CHECK(text.contains("> <CHARGEFW_CHARGES_2>\n-0.7000 0.3500 0.3500\n\n"));
        CHECK(text.ends_with("$$$$\n"));
    }
}
