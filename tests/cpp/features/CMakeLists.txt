chargefw_add_snitch_test(test_features
        test_topology_features.cpp
        test_conformer_features.cpp
        test_prepared_molecule.cpp
        test_spatial_fragment.cpp
)

target_include_directories(test_features PRIVATE
        ${PROJECT_SOURCE_DIR}/src
)
