chargefw_add_snitch_test(test_topology_features
        test_topology_features.cpp
)

target_include_directories(test_topology_features PRIVATE
        ${PROJECT_SOURCE_DIR}/src
)

chargefw_add_snitch_test(test_conformer_features
        test_conformer_features.cpp
)

chargefw_add_snitch_test(test_prepared_molecule
        test_prepared_molecule.cpp
)

chargefw_add_snitch_test(test_spatial_fragment
        test_spatial_fragment.cpp
)
