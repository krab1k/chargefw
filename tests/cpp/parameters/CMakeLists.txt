chargefw_add_snitch_test(test_parameter_classifier
        test_parameter_classifier.cpp
)

chargefw_add_snitch_test(test_parameter_set_io
        test_parameter_set_io.cpp
)

target_compile_definitions(test_parameter_set_io
        PRIVATE
        CHARGEFW_TEST_PARAMETER_DIR="${PROJECT_SOURCE_DIR}/data/parameters"
)