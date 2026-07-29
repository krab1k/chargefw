execute_process(
        COMMAND "${CMAKE_COMMAND}" -E env
                "CHARGEFW_PARAMETER_DIR=${CHARGEFW_PARAMETER_DIR}"
                "${CHARGEFW_CLI}" "${CHARGEFW_INPUT}"
        RESULT_VARIABLE result
        OUTPUT_VARIABLE output
        ERROR_VARIABLE error
)

if(NOT result EQUAL 0)
    message(FATAL_ERROR "chargefw failed with exit code ${result}: ${error}")
endif()

string(JSON schema_version GET "${output}" schema_version)
if(NOT schema_version STREQUAL "1.0")
    message(FATAL_ERROR "Unexpected schema version: ${schema_version}")
endif()

string(JSON generator_name GET "${output}" generator name)
if(NOT generator_name STREQUAL "ChargeFW")
    message(FATAL_ERROR "Unexpected generator name: ${generator_name}")
endif()

string(JSON result_count LENGTH "${output}" results)
if(NOT result_count EQUAL 1)
    message(FATAL_ERROR "Expected one result, got ${result_count}")
endif()

string(JSON status GET "${output}" results 0 status)
if(NOT status STREQUAL "success")
    message(FATAL_ERROR "Expected successful calculation, got ${status}")
endif()

string(JSON atom_mapping_kind GET "${output}" results 0 input atom_mapping kind)
if(NOT atom_mapping_kind STREQUAL "identity")
    message(FATAL_ERROR "Expected identity atom mapping, got ${atom_mapping_kind}")
endif()

string(JSON conformer_mapping_kind GET "${output}" results 0 input conformer_mapping kind)
if(NOT conformer_mapping_kind STREQUAL "identity")
    message(FATAL_ERROR "Expected identity conformer mapping, got ${conformer_mapping_kind}")
endif()

string(JSON assignment_count LENGTH "${output}" results 0 assignments)
if(NOT assignment_count EQUAL 1)
    message(FATAL_ERROR "Expected one assignment, got ${assignment_count}")
endif()

string(JSON atom_order GET "${output}" results 0 assignments 0 atom_order)
if(NOT atom_order STREQUAL "source")
    message(FATAL_ERROR "Expected source atom order, got ${atom_order}")
endif()

string(JSON charge_count LENGTH "${output}" results 0 assignments 0 charges)
if(NOT charge_count EQUAL 3)
    message(FATAL_ERROR "Expected three charges, got ${charge_count}")
endif()
