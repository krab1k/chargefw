#include <chargefw/core/atom.h>
#include <chargefw/core/bond.h>
#include <chargefw/core/conformer.h>
#include <chargefw/core/molecule.h>
#include <chargefw/core/position.h>

#include <chargefw/features/conformer_features.h>
#include <chargefw/features/topology_features.h>

#include <chargefw/methods/calculation_input.h>
#include <chargefw/methods/method_options.h>
#include <chargefw/methods/method_registry.h>

#include <chargefw/parameters/atom_parameters.h>
#include <chargefw/parameters/common_parameters.h>
#include <chargefw/parameters/parameter_classifier.h>
#include <chargefw/parameters/parameter_key.h>
#include <chargefw/parameters/parameter_set.h>
#include <chargefw/parameters/parameter_view.h>

#include <iostream>
#include <stdexcept>
#include <string>
#include <string_view>
#include <utility>
#include <vector>

namespace core = chargefw::core;
namespace features = chargefw::features;
namespace methods = chargefw::methods;
namespace params = chargefw::parameters;

auto make_water() -> core::Molecule
{
    auto atoms = std::vector{
        core::Atom{8, 0, "O"},
        core::Atom{1, 0, "H1"},
        core::Atom{1, 0, "H2"}
    };

    auto bonds = std::vector{
        core::Bond{0, 1, core::BondOrder::SINGLE},
        core::Bond{0, 2, core::BondOrder::SINGLE}
    };

    auto positions = std::vector{
        core::Position{0.0000, 0.0000, 0.0000},
        core::Position{0.9572, 0.0000, 0.0000},
        core::Position{-0.2390, 0.9270, 0.0000}
    };

    auto conformers = std::vector{
        core::Conformer{std::move(positions), "model-1"}
    };

    return core::Molecule{
        std::move(atoms),
        std::move(bonds),
        std::move(conformers),
        "water"
    };
}

auto make_example_eem_parameters() -> params::ParameterSet
{
    return params::ParameterSet{
        params::ParameterSetMetadata{
            .id = "example-eem",
            .method_id = "eem",
            .name = "Example EEM",
            .publication = "internal",
            .notes = "development fixture"
        },
        params::CommonParameters{
            std::vector{
                params::NamedParameter{.name = "kappa", .value = 1.0}
            }
        },
        params::AtomParameters{
            std::vector{
                params::AtomParameterEntry{
                    .key = params::AtomParameterKey{
                        .atomic_number = 1,
                        .classification = params::AtomParameterClassificationKind::PLAIN,
                        .type = "*"
                    },
                    .parameters = std::vector{
                        params::NamedParameter{.name = "A", .value = 7.17},
                        params::NamedParameter{.name = "B", .value = 13.89}
                    }
                },
                params::AtomParameterEntry{
                    .key = params::AtomParameterKey{
                        .atomic_number = 8,
                        .classification = params::AtomParameterClassificationKind::PLAIN,
                        .type = "*"
                    },
                    .parameters = std::vector{
                        params::NamedParameter{.name = "A", .value = 8.74},
                        params::NamedParameter{.name = "B", .value = 13.36}
                    }
                }
            }
        },
        params::BondParameters{}
    };
}

auto require_method(const methods::MethodRegistry& registry, std::string_view id)
    -> const methods::Method&
{
    const auto* method = registry.find(id);

    if (method == nullptr) {
        throw std::runtime_error{"method not found: " + std::string{id}};
    }

    return *method;
}

auto calculate(const methods::Method& method,
               const core::Molecule& molecule,
               const features::TopologyFeatures& topology,
               const features::ConformerFeatures* geometry)
    -> chargefw::charges::AtomicCharges
{
    const auto method_options = methods::make_default_options(method.option_schema());

    const auto input = methods::CalculationInput{
        .molecule = molecule,
        .topology = topology,
        .geometry = geometry,
        .method_options = method_options
    };

    return method.calculate(input);
}

auto print_charges(std::string_view title,
                   const core::Molecule& molecule,
                   const chargefw::charges::AtomicCharges& charges) -> void
{
    std::cout << title << '\n';

    for (std::size_t atom_index = 0; atom_index < molecule.atom_count(); ++atom_index) {
        const auto& atom = molecule.atom(atom_index);

        std::cout
            << "  atom " << atom_index
            << " " << atom.name()
            << " Z=" << atom.atomic_number()
            << " q=" << charges[atom_index]
            << '\n';
    }

    std::cout << "  total = " << charges.total() << "\n\n";
}

auto main() -> int
{
    const auto molecule = make_water();

    const features::TopologyFeatures topology{molecule};
    const features::ConformerFeatures geometry{molecule, 0};

    std::cout << "Molecule: " << molecule.name() << '\n';
    std::cout << "Atoms: " << molecule.atom_count() << '\n';
    std::cout << "Bonds: " << molecule.bond_count() << "\n\n";

    std::cout << "Topology features:\n";
    for (std::size_t atom_index = 0; atom_index < molecule.atom_count(); ++atom_index) {
        std::cout
            << "  atom " << atom_index
            << " degree=" << topology.degree(atom_index)
            << " neighbors=";

        for (const auto neighbor_index : topology.neighbor_indices(atom_index)) {
            std::cout << neighbor_index << ' ';
        }

        std::cout << '\n';
    }

    std::cout << "\nGeometry features:\n";
    std::cout << "  O-H1 distance = " << geometry.distance(0, 1) << '\n';
    std::cout << "  O-H2 distance = " << geometry.distance(0, 2) << "\n\n";

    const auto& registry = methods::method_registry();

    const auto& dummy = require_method(registry, "dummy");
    const auto& formal = require_method(registry, "formal");

    const auto dummy_charges = calculate(dummy, molecule, topology, &geometry);
    const auto formal_charges = calculate(formal, molecule, topology, &geometry);

    print_charges("Dummy charges:", molecule, dummy_charges);
    print_charges("Formal charges:", molecule, formal_charges);

    const auto parameter_set = make_example_eem_parameters();

    const auto classification = params::classify_parameters(
        molecule,
        topology,
        parameter_set
    );

    const auto parameters = params::ParameterView{
        parameter_set,
        classification
    };

    const auto atom_a = parameters.atom("A");
    const auto atom_b = parameters.atom("B");
    const auto kappa = parameters.common("kappa");

    std::cout << "Parameter view:\n";
    std::cout << "  parameter set = " << parameter_set.id() << '\n';
    std::cout << "  kappa = " << kappa << '\n';

    for (std::size_t atom_index = 0; atom_index < molecule.atom_count(); ++atom_index) {
        const auto& atom = molecule.atom(atom_index);

        std::cout
            << "  atom " << atom_index
            << " " << atom.name()
            << " A=" << atom_a[atom_index]
            << " B=" << atom_b[atom_index]
            << '\n';
    }

    return 0;
}