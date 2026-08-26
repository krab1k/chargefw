#include <chargefw/adapters/native/mol2_output.h>
#include <chargefw/charges/atomic_charges.h>
#include <chargefw/charges/charge_collection.h>
#include <chargefw/core/atom.h>
#include <chargefw/core/bond.h>
#include <chargefw/core/conformer.h>
#include <chargefw/core/molecule.h>
#include <chargefw/core/position.h>
#include <snitch/snitch.hpp>

#include <filesystem>
#include <fstream>
#include <print>
#include <sstream>
#include <string>
#include <vector>

namespace mol2_output = chargefw::adapters::native::mol2_output;
namespace charges = chargefw::charges;

namespace {

[[nodiscard]] auto fixture(const std::string_view name) -> std::filesystem::path {
    return std::filesystem::path{CHARGEFW_TEST_SOURCE_DIR} / "tests" / "fixtures" / name;
}

[[nodiscard]] auto assignment(std::vector<double> values) -> charges::ChargeAssignment {
    return {.target = {.molecule_index = 0, .conformer_index = 0},
            .charges = charges::AtomicCharges{std::move(values)}};
}

[[nodiscard]] auto assignment(const std::size_t molecule_index, std::vector<double> values)
    -> charges::ChargeAssignment {
    return {.target = {.molecule_index = molecule_index, .conformer_index = 0},
            .charges = charges::AtomicCharges{std::move(values)}};
}

[[nodiscard]] auto geometry_independent_assignment(std::vector<double> values)
    -> charges::ChargeAssignment {
    return {.target = {.molecule_index = 0}, .charges = charges::AtomicCharges{std::move(values)}};
}

} // namespace

TEST_CASE("MOL2 output preserves source structure and replaces atom charges", "[adapters][mol2]") {
    {
        auto output = std::ostringstream{};
        auto input = std::ifstream{fixture("synthetic/mol2/aromatic.mol2"), std::ios::binary};
        const auto source = std::string{std::istreambuf_iterator<char>{input}, {}};
        const auto assignments = std::vector{assignment({-0.87654, 0.43827, 0.43827})};
        mol2_output::Mol2Writer{output}.write_preserving_buffer(source, assignments);
        const auto text = output.str();
        CHECK(text.contains("7 N1 0.0000 0.0000 0.0000 N.4 1 LIG -0.8765"));
        CHECK(text.contains("9 C1 1.0000 0.0000 0.0000 C.ar 1 LIG 0.4383"));
        CHECK(text.contains("11 O1 2.0000 0.0000 0.0000 O.2 1 LIG 0.4383"));
        CHECK(text.contains("@<TRIPOS>BOND\n1 7 9 ar\n2 9 11 2\n"));

        auto expected = source;
        const auto first_charge = expected.find("0.2500", expected.find("1 LIG"));
        expected.replace(first_charge, std::string_view{"0.2500"}.size(), "-0.8765");
        const auto second_charge =
            expected.find("0.0000", expected.find("1 LIG", first_charge + 1));
        expected.replace(second_charge, std::string_view{"0.0000"}.size(), "0.4383");
        const auto third_charge =
            expected.find("0.0000", expected.find("1 LIG", second_charge + 1));
        expected.replace(third_charge, std::string_view{"0.0000"}.size(), "0.4383");
        CHECK(text == expected);
    }

    {
        const auto source =
            std::filesystem::path{CHARGEFW_TEST_SOURCE_DIR} / "build" / "mol2_without_charges.mol2";
        auto input = std::ofstream{source};
        std::print(input, "@<TRIPOS>MOLECULE\nminimal\n2 0 0 0 0\nSMALL\n"
                          "NO_CHARGES\n\n@<TRIPOS>ATOM\n1 C1 0 0 0 C.3\n"
                          "2 H1 1 0 0 H\n@<TRIPOS>BOND\n");
        input.close();

        auto output = std::ostringstream{};
        const auto assignments = std::vector{assignment({-0.1, 0.1})};
        mol2_output::Mol2Writer{output}.write_preserving_source(source.string(), assignments);
        const auto text = output.str();
        CHECK(text.contains("1 C1 0 0 0 C.3 1 UNL -0.1000"));
        CHECK(text.contains("2 H1 1 0 0 H 1 UNL 0.1000"));
        std::filesystem::remove(source);
    }

    {
        const auto source = std::filesystem::path{CHARGEFW_TEST_SOURCE_DIR} / "build" /
                            "mol2_batch_preserving.mol2";
        auto input = std::ofstream{source};
        std::print(input, "preamble\n@<TRIPOS>MOLECULE\nfirst\n1 0 0 0 0\nSMALL\nNO_CHARGES\n"
                          "@<TRIPOS>ATOM\n1 C1 0 0 0 C.3\n@<TRIPOS>BOND\n"
                          "@<TRIPOS>MOLECULE\nsecond\n1 0 0 0 0\nSMALL\nUSER_CHARGES\n"
                          "@<TRIPOS>ATOM\n1 O1 0 0 0 O.3 1 UNL 0.0000\n@<TRIPOS>BOND\n");
        input.close();

        auto output = std::ostringstream{};
        const auto assignments = std::vector{assignment(0, {-0.1}), assignment(1, {0.2})};
        mol2_output::Mol2Writer{output}.write_preserving_source(source.string(), assignments);
        CHECK(output.str().contains("1 C1 0 0 0 C.3 1 UNL -0.1000"));
        CHECK(output.str().contains("1 O1 0 0 0 O.3 1 UNL 0.2000"));
        CHECK(output.str().contains("preamble\n"));
        std::filesystem::remove(source);
    }

    {
        const auto source = std::filesystem::path{CHARGEFW_TEST_SOURCE_DIR} / "build" /
                            "mol2_preserved_spacing.mol2";
        auto input = std::ofstream{source, std::ios::binary};
        std::print(input, "@<TRIPOS>MOLECULE\r\nminimal\r\n1 0 0 0 0\r\nSMALL\r\n"
                          "USER_CHARGES\r\n\r\n@<TRIPOS>ATOM\r\n"
                          "  1\tC1   0.0  0.0\t0.0  C.3  1 LIG   0.1234   # retained\r\n"
                          "@<TRIPOS>BOND\r\n");
        input.close();

        auto output = std::ostringstream{};
        const auto assignments = std::vector{assignment({-0.1})};
        mol2_output::Mol2Writer{output}.write_preserving_source(source.string(), assignments);
        CHECK(output.str() == "@<TRIPOS>MOLECULE\r\nminimal\r\n1 0 0 0 0\r\nSMALL\r\n"
                              "USER_CHARGES\r\n\r\n@<TRIPOS>ATOM\r\n"
                              "  1\tC1   0.0  0.0\t0.0  C.3  1 LIG   -0.1000   # retained\r\n"
                              "@<TRIPOS>BOND\r\n");
        std::filesystem::remove(source);
    }

    {
        const auto molecule = chargefw::core::Molecule{
            {chargefw::core::Atom{8}, chargefw::core::Atom{1}, chargefw::core::Atom{1}},
            {chargefw::core::Bond{0, 1}, chargefw::core::Bond{0, 2}},
            {chargefw::core::Conformer{{{.x = 0.0, .y = 0.0, .z = 0.0},
                                        {.x = 0.9, .y = 0.0, .z = 0.0},
                                        {.x = -0.2, .y = 0.9, .z = 0.0}}}},
            "water"};
        auto output = std::ostringstream{};
        mol2_output::Mol2Writer{output}.write_generated(molecule,
                                                        assignment({-0.97533, 0.48766, 0.48767}));
        const auto text = output.str();
        CHECK(text.contains("@<TRIPOS>MOLECULE\nwater\n3 2 0 0 0\nSMALL\nUSER_CHARGES\n"));
        CHECK(text.contains("1 O1 0 0 0 O 1 UNL -0.9753"));
        CHECK(text.contains("2 H2 0.9 0 0 H 1 UNL 0.4877"));
        CHECK(text.contains("@<TRIPOS>BOND\n1 1 2 1\n2 1 3 1\n"));
    }

    {
        const auto molecule = chargefw::core::Molecule{
            {chargefw::core::Atom{8}, chargefw::core::Atom{1}, chargefw::core::Atom{1}},
            {chargefw::core::Bond{0, 1}, chargefw::core::Bond{0, 2}},
            {chargefw::core::Conformer{{{.x = 0.0, .y = 0.0, .z = 0.0},
                                        {.x = 0.9, .y = 0.0, .z = 0.0},
                                        {.x = -0.2, .y = 0.9, .z = 0.0}}}},
            "water"};
        auto output = std::ostringstream{};
        mol2_output::Mol2Writer{output}.write_generated(
            molecule, geometry_independent_assignment({-0.97533, 0.48766, 0.48767}));
        CHECK(output.str().contains("1 O1 0 0 0 O 1 UNL -0.9753"));
    }
}
