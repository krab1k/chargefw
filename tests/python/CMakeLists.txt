add_test(
    NAME test_chargefw_python_import
    COMMAND "${CHARGEFW_RESOLVED_PYTHON_EXECUTABLE}" "${CMAKE_CURRENT_SOURCE_DIR}/test_import.py"
            "${PROJECT_VERSION}"
)
set_tests_properties(test_chargefw_python_import PROPERTIES
    ENVIRONMENT "PYTHONPATH=${PROJECT_BINARY_DIR}/python"
)

add_test(
    NAME test_chargefw_python_molecule
    COMMAND "${CHARGEFW_RESOLVED_PYTHON_EXECUTABLE}" "${CMAKE_CURRENT_SOURCE_DIR}/test_molecule.py"
)
set_tests_properties(test_chargefw_python_molecule PROPERTIES
    ENVIRONMENT "PYTHONPATH=${PROJECT_BINARY_DIR}/python"
)

add_test(
    NAME test_chargefw_python_calculation
    COMMAND "${CHARGEFW_RESOLVED_PYTHON_EXECUTABLE}" "${CMAKE_CURRENT_SOURCE_DIR}/test_calculation.py"
)
set_tests_properties(test_chargefw_python_calculation PROPERTIES
    ENVIRONMENT "PYTHONPATH=${PROJECT_BINARY_DIR}/python"
)
