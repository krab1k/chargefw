chargefw_add_test(test_builtin_methods
        test_builtin_methods.cpp
)

chargefw_add_test(test_method_options
        test_method_options.cpp
)

chargefw_add_test(test_method_registry
        test_method_registry.cpp
)

chargefw_add_test(test_method_prerequisites
        test_method_prerequisites.cpp
)

chargefw_add_test(test_method_applicability
        test_method_applicability.cpp
)

chargefw_add_test(test_method_calculation
        test_method_calculation.cpp
)

chargefw_add_test(test_peoe
        test_peoe.cpp
)

chargefw_add_test(test_mpeoe
        test_mpeoe.cpp
)

target_compile_definitions(test_mpeoe
        PRIVATE
        CHARGEFW_TEST_PARAMETER_DIR="${PROJECT_SOURCE_DIR}/data/parameters"
)

chargefw_add_test(test_veem
        test_veem.cpp
)

chargefw_add_test(test_gdac
        test_gdac.cpp
)

target_compile_definitions(test_gdac
        PRIVATE
        CHARGEFW_TEST_PARAMETER_DIR="${PROJECT_SOURCE_DIR}/data/parameters"
)

chargefw_add_test(test_delre
        test_delre.cpp
)

chargefw_add_test(test_denr
        test_denr.cpp
)

chargefw_add_test(test_kcm
        test_kcm.cpp
)

chargefw_add_test(test_qeq
        test_qeq.cpp
)

chargefw_add_test(test_eem
        test_eem.cpp
)

chargefw_add_test(test_smpqeq
        test_smpqeq.cpp
)

chargefw_add_test(test_sfkeem
        test_sfkeem.cpp
)