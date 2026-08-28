add_test(
    NAME test_chargefw_python_import
    COMMAND "${CHARGEFW_RESOLVED_PYTHON_EXECUTABLE}" "${CMAKE_CURRENT_SOURCE_DIR}/test_import.py"
)
set_tests_properties(test_chargefw_python_import PROPERTIES
    ENVIRONMENT "PYTHONPATH=${PROJECT_BINARY_DIR}/python;CHARGEFW_EXPECTED_VERSION=${PROJECT_VERSION}"
)

add_test(
    NAME test_chargefw_python_component_install
    COMMAND ${CMAKE_COMMAND}
            -DCHARGEFW_BINARY_DIR=${PROJECT_BINARY_DIR}
            -DCHARGEFW_EXPECTED_VERSION=${PROJECT_VERSION}
            -DCHARGEFW_PYTHON_EXECUTABLE=${CHARGEFW_RESOLVED_PYTHON_EXECUTABLE}
            -DCHARGEFW_PYTHON_INSTALL_PACKAGE_DIRECTORY=${CHARGEFW_PYTHON_INSTALL_PACKAGE_DIRECTORY}
            -DCHARGEFW_TEST_SOURCE_DIR=${CMAKE_CURRENT_SOURCE_DIR}
            -P ${CMAKE_CURRENT_SOURCE_DIR}/test_chargefw_python_component_install.cmake
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

add_test(
    NAME test_chargefw_python_mypy
    COMMAND "${CHARGEFW_RESOLVED_PYTHON_EXECUTABLE}" -m mypy --strict python/chargefw tests/python
)
set_tests_properties(test_chargefw_python_mypy PROPERTIES
    WORKING_DIRECTORY "${PROJECT_SOURCE_DIR}"
    ENVIRONMENT "MYPYPATH=${PROJECT_SOURCE_DIR}/python"
)
