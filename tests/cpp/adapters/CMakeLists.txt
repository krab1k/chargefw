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

chargefw_add_snitch_test(test_pdb
        test_pdb.cpp
)

chargefw_add_snitch_test(test_mmcif
        test_mmcif.cpp
)

chargefw_add_snitch_test(test_mmcif_output
        test_mmcif_output.cpp
)

chargefw_add_snitch_test(test_fixed_ions
        test_fixed_ions.cpp
)

target_include_directories(test_fixed_ions PRIVATE ${PROJECT_SOURCE_DIR}/src)
target_compile_definitions(test_fixed_ions PRIVATE
        CHARGEFW_TEST_PARAMETER_DIR="${PROJECT_SOURCE_DIR}/data/parameters"
)

target_compile_definitions(test_mol
        PRIVATE
        CHARGEFW_TEST_SOURCE_DIR="${PROJECT_SOURCE_DIR}"
)
