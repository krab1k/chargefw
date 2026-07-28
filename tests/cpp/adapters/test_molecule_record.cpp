#include "support/test_molecules.h"

#include <cassert>
#include <chargefw/adapters/molecule_record.h>
#include <optional>

namespace adapters = chargefw::adapters;

auto main() -> int {
    const auto water = chargefw::test::make_water();

    const adapters::MoleculeRecordMapping identity{.atom_indices = {0, 1, 2},
                                                   .conformer_indices = {0}};
    assert(adapters::is_identity_mapping(identity, water));

    const adapters::MoleculeRecordMapping reordered{.atom_indices = {0, 2, 1},
                                                    .conformer_indices = {0}};
    assert(!adapters::is_identity_mapping(reordered, water));

    const adapters::MoleculeRecordMapping omitted{.atom_indices = {0, std::nullopt, 2},
                                                  .conformer_indices = {0}};
    assert(!adapters::is_identity_mapping(omitted, water));

    const adapters::MoleculeRecordMapping incomplete{.atom_indices = {0, 1, 2},
                                                     .conformer_indices = {}};
    assert(!adapters::is_identity_mapping(incomplete, water));

    const adapters::ImportedMoleculeRecord imported{
        .molecule = chargefw::test::make_water(),
        .identity = {.source = "input.cif", .record_index = 0, .record_id = "example"},
        .mapping = identity,
        .diagnostics = {},
        .source = adapters::MoleculeRecordSource{.format = "mmcif", .payload = "data_example\n"}};
    assert(imported.source.has_value());
    assert(imported.source->format == "mmcif");
    assert(imported.source->payload == "data_example\n");

    return 0;
}
