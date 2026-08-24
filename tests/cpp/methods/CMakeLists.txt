chargefw_add_snitch_test(test_builtin_methods
        test_builtin_methods.cpp
)

target_compile_definitions(test_builtin_methods
        PRIVATE
        CHARGEFW_TEST_PARAMETER_DIR="${PROJECT_SOURCE_DIR}/data/parameters"
)

chargefw_add_snitch_test(test_method_options
        test_method_options.cpp
)

chargefw_add_snitch_test(test_method_registry
        test_method_registry.cpp
)

chargefw_add_snitch_test(test_method_prerequisites
        test_method_prerequisites.cpp
)

chargefw_add_snitch_test(test_method_applicability
        test_method_applicability.cpp
)

chargefw_add_snitch_test(test_method_calculation
        test_method_calculation.cpp
)

chargefw_add_snitch_test(test_calculation_input
        test_calculation_input.cpp
)

chargefw_add_snitch_test(test_peoe
        test_peoe.cpp
)

chargefw_add_snitch_test(test_mpeoe
        test_mpeoe.cpp
)

chargefw_add_snitch_test(test_veem
        test_veem.cpp
)

chargefw_add_snitch_test(test_gdac
        test_gdac.cpp
)

chargefw_add_snitch_test(test_tsef
        test_tsef.cpp
)

chargefw_add_snitch_test(test_delre
        test_delre.cpp
)

chargefw_add_snitch_test(test_denr
        test_denr.cpp
)

chargefw_add_snitch_test(test_kcm
        test_kcm.cpp
)

chargefw_add_snitch_test(test_qeq
        test_qeq.cpp
)

chargefw_add_snitch_test(test_sqe
        test_sqe.cpp
)

chargefw_add_snitch_test(test_sqeq0
        test_sqeq0.cpp
)

chargefw_add_snitch_test(test_sqeqp
        test_sqeqp.cpp
)

chargefw_add_snitch_test(test_eem
        test_eem.cpp
)

chargefw_add_snitch_test(test_smpqeq
        test_smpqeq.cpp
)

chargefw_add_snitch_test(test_sfkeem
        test_sfkeem.cpp
)

chargefw_add_snitch_test(test_eqeq
        test_eqeq.cpp
)

chargefw_add_snitch_test(test_eqeqc
        test_eqeqc.cpp
)

chargefw_add_snitch_test(test_abeem
        test_abeem.cpp
)
