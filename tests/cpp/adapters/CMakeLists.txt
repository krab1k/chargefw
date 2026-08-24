chargefw_add_snitch_test(test_mol
        test_mol.cpp
)

chargefw_add_snitch_test(test_json
        test_json.cpp
)

chargefw_add_snitch_test(test_json_output
        test_json_output.cpp
)

chargefw_add_snitch_test(test_mol2_output
        test_mol2_output.cpp
)

chargefw_add_snitch_test(test_sdf_output
        test_sdf_output.cpp
)

chargefw_add_snitch_test(test_pdb
        test_pdb.cpp
)

chargefw_add_snitch_test(test_mmcif
        test_mmcif.cpp
)

chargefw_add_snitch_test(test_mmcif_output
        test_mmcif_output.cpp
)

target_compile_definitions(test_mol
        PRIVATE
        CHARGEFW_TEST_SOURCE_DIR="${PROJECT_SOURCE_DIR}"
)

target_compile_definitions(test_mol2_output
        PRIVATE
        CHARGEFW_TEST_SOURCE_DIR="${PROJECT_SOURCE_DIR}"
)
