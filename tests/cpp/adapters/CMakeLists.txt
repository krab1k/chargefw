chargefw_add_test(test_molecule_record
        test_molecule_record.cpp
)

chargefw_add_test(test_mol
        test_mol.cpp
)

chargefw_add_test(test_json
        test_json.cpp
)

chargefw_add_test(test_json_output
        test_json_output.cpp
)

chargefw_add_test(test_mol2_output
        test_mol2_output.cpp
)

target_compile_definitions(test_mol
        PRIVATE
        CHARGEFW_TEST_SOURCE_DIR="${PROJECT_SOURCE_DIR}"
)

target_compile_definitions(test_mol2_output
        PRIVATE
        CHARGEFW_TEST_SOURCE_DIR="${PROJECT_SOURCE_DIR}"
)
