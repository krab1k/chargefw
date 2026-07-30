#include <cassert>
#include <chargefw/adapters/native/mol2_output.h>
#include <chargefw/charges/atomic_charges.h>
#include <chargefw/charges/charge_collection.h>
#include <chargefw/core/atom.h>
#include <chargefw/core/bond.h>
#include <chargefw/core/conformer.h>
#include <chargefw/core/molecule.h>
#include <chargefw/core/position.h>

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
    return std::filesystem::path{CHARGEFW_TEST_SOURCE_DIR} / "tests" / "fixtures" / "mol2" / name;
}

[[nodiscard]] auto assignment(std::vector<double> values) -> charges::ChargeAssignment {
    return {.target = {.molecule_index = 0, .conformer_index = 0},
            .charges = charges::AtomicCharges{std::move(values)}};
}

} // namespace

auto main() -> int {
    {
        auto output = std::ostringstream{};
        mol2_output::Mol2Writer{output}.write_preserving_source(
            fixture("charged_aromatic.mol2").string(), 0, assignment({-0.87654, 0.43827, 0.43827}));
        const auto text = output.str();
        assert(text.contains("7 N1 0.0000 0.0000 0.0000 N.4 1 LIG -0.8765"));
        assert(text.contains("9 C1 1.0000 0.0000 0.0000 C.ar 1 LIG 0.4383"));
        assert(text.contains("11 O1 2.0000 0.0000 0.0000 O.2 1 LIG 0.4383"));
        assert(text.contains("@<TRIPOS>BOND\n1 7 9 ar\n2 9 11 2\n"));
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
        mol2_output::Mol2Writer{output}.write_preserving_source(source.string(), 0,
                                                                assignment({-0.1, 0.1}));
        const auto text = output.str();
        assert(text.contains("1 C1 0 0 0 C.3 1 UNL -0.1000"));
        assert(text.contains("2 H1 1 0 0 H 1 UNL 0.1000"));
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
        mol2_output::Mol2Writer{output}.write_preserving_source(source.string(), 0,
                                                                assignment({-0.1}));
        assert(output.str() == "@<TRIPOS>MOLECULE\r\nminimal\r\n1 0 0 0 0\r\nSMALL\r\n"
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
        assert(text.contains("@<TRIPOS>MOLECULE\nwater\n3 2 0 0 0\nSMALL\nUSER_CHARGES\n"));
        assert(text.contains("1 O1 0 0 0 O 1 UNL -0.9753"));
        assert(text.contains("2 H2 0.9 0 0 H 1 UNL 0.4877"));
        assert(text.contains("@<TRIPOS>BOND\n1 1 2 1\n2 1 3 1\n"));
    }

    return 0;
}
