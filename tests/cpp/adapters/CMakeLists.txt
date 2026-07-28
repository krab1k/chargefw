chargefw_add_test(test_molecule_record
        test_molecule_record.cpp
)

chargefw_add_test(test_mol
        test_mol.cpp
)

target_compile_definitions(test_mol
        PRIVATE
        CHARGEFW_TEST_SOURCE_DIR="${PROJECT_SOURCE_DIR}"
)
