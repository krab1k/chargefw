#include <chargefw/adapters/native/mol.h>
#include <chargefw/adapters/native/sdf.h>

#include <istream>
#include <string>
#include <utility>

namespace chargefw::adapters::native::sdf {
namespace {

auto consume_to_sdf_delimiter(std::istream& input) -> void {
    std::string line;

    while (std::getline(input, line)) {
        if (line == "$$$$") {
            return;
        }
    }
}

} // namespace

SdfReader::SdfReader(std::istream& input, std::string source)
    : input_{std::addressof(input)}, source_{std::move(source)} {}

auto SdfReader::next() -> std::optional<::chargefw::adapters::MoleculeRecordResult> {
    if (input_->peek() == std::char_traits<char>::eof()) {
        return std::nullopt;
    }

    const auto identity =
        MoleculeRecordIdentity{.source = source_, .record_index = record_index_++, .record_id = {}};
    auto result = mol::parse_mol(*input_, identity);
    consume_to_sdf_delimiter(*input_);
    return result;
}

} // namespace chargefw::adapters::native::sdf
