set(CHARGEFW_INSTALL_PREFIX "${CMAKE_CURRENT_BINARY_DIR}/chargefw_python_component_prefix")
file(REMOVE_RECURSE "${CHARGEFW_INSTALL_PREFIX}")

execute_process(
    COMMAND "${CMAKE_COMMAND}" --install "${CHARGEFW_BINARY_DIR}" --component python
            --prefix "${CHARGEFW_INSTALL_PREFIX}"
    RESULT_VARIABLE install_result
    ERROR_VARIABLE install_error
)
if(NOT install_result EQUAL 0)
    message(FATAL_ERROR "ChargeFW Python component installation failed: ${install_error}")
endif()

get_filename_component(CHARGEFW_PYTHON_SITEARCH
                       "${CHARGEFW_PYTHON_INSTALL_PACKAGE_DIRECTORY}" DIRECTORY)
if(EXISTS "${CHARGEFW_INSTALL_PREFIX}/chargefw")
    message(FATAL_ERROR "ChargeFW Python package was installed outside site-packages")
endif()
if(NOT EXISTS "${CHARGEFW_INSTALL_PREFIX}/${CHARGEFW_PYTHON_INSTALL_PACKAGE_DIRECTORY}/__init__.py")
    message(FATAL_ERROR "ChargeFW Python package was not installed into site-packages")
endif()

execute_process(
    COMMAND "${CMAKE_COMMAND}" -E env "PYTHONPATH=${CHARGEFW_INSTALL_PREFIX}/${CHARGEFW_PYTHON_SITEARCH}"
            "CHARGEFW_EXPECTED_VERSION=${CHARGEFW_EXPECTED_VERSION}"
            "${CHARGEFW_PYTHON_EXECUTABLE}" "${CHARGEFW_TEST_SOURCE_DIR}/test_import.py"
    RESULT_VARIABLE import_result
    ERROR_VARIABLE import_error
)
if(NOT import_result EQUAL 0)
    message(FATAL_ERROR "Installed ChargeFW Python package import failed: ${import_error}")
endif()
