set(output_directory "${CMAKE_CURRENT_BINARY_DIR}/chargefw_cli_outputs")
set(output_prefix "${output_directory}/water.chargefw")
file(REMOVE_RECURSE "${output_directory}")

execute_process(
        COMMAND "${CMAKE_COMMAND}" -E env
                "CHARGEFW_PARAMETER_DIR=${CHARGEFW_PARAMETER_DIR}"
                "${CHARGEFW_CLI}" "${CHARGEFW_INPUT}" "${output_directory}"
        RESULT_VARIABLE result
        OUTPUT_VARIABLE output
        ERROR_VARIABLE error
)

if(NOT result EQUAL 0)
    message(FATAL_ERROR "chargefw failed with exit code ${result}: ${error}")
endif()

if(NOT IS_DIRECTORY "${output_directory}")
    message(FATAL_ERROR "Output directory was not created: ${output_directory}")
endif()

foreach(extension IN ITEMS json sdf mol2 cif)
    set(path "${output_prefix}.${extension}")
    if(NOT EXISTS "${path}")
        message(FATAL_ERROR "Expected output file was not created: ${path}")
    endif()
endforeach()

file(READ "${output_prefix}.json" json_output)
string(JSON schema_version GET "${json_output}" schema_version)
if(NOT schema_version STREQUAL "1.0")
    message(FATAL_ERROR "Unexpected schema version: ${schema_version}")
endif()
string(JSON status GET "${json_output}" results 0 status)
if(NOT status STREQUAL "success")
    message(FATAL_ERROR "Expected successful calculation, got ${status}")
endif()

file(READ "${output_prefix}.sdf" sdf_output)
if(NOT sdf_output MATCHES "> <CHARGEFW_CHARGES_1>")
    message(FATAL_ERROR "SDF output does not contain ChargeFW charges")
endif()
if(NOT sdf_output MATCHES "PUBCHEM_COMPOUND_CID")
    message(FATAL_ERROR "SDF source properties were not preserved")
endif()

file(READ "${output_prefix}.mol2" mol2_output)
if(NOT mol2_output MATCHES "@<TRIPOS>MOLECULE")
    message(FATAL_ERROR "MOL2 output was not generated")
endif()
if(NOT mol2_output MATCHES "USER_CHARGES")
    message(FATAL_ERROR "MOL2 output does not declare generated charges")
endif()

file(REMOVE_RECURSE "${output_directory}")
