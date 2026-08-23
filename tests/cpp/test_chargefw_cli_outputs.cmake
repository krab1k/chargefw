set(output_directory "${CMAKE_CURRENT_BINARY_DIR}/chargefw_cli_outputs")
set(output_prefix "${output_directory}/water.chargefw")
file(REMOVE_RECURSE "${output_directory}")

foreach(mode IN ITEMS full cutoff cover)
    set(mode_output_directory "${CMAKE_CURRENT_BINARY_DIR}/chargefw_cli_${mode}")
    set(mode_output_prefix "${mode_output_directory}/water.chargefw")
    file(REMOVE_RECURSE "${mode_output_directory}")

    set(execution_arguments --method eem --execution ${mode})
    if(NOT mode STREQUAL "full")
        list(APPEND execution_arguments --radius 8)
    endif()
    if(mode STREQUAL "cutoff")
        list(APPEND execution_arguments --progress)
    endif()

    execute_process(
            COMMAND "${CMAKE_COMMAND}" -E env
                    "CHARGEFW_PARAMETER_DIR=${CHARGEFW_PARAMETER_DIR}"
                    "${CHARGEFW_CLI}" calculate ${execution_arguments} "${CHARGEFW_INPUT}"
                    "${mode_output_directory}"
            RESULT_VARIABLE mode_result
            ERROR_VARIABLE mode_error
    )
    if(NOT mode_result EQUAL 0)
        message(FATAL_ERROR "${mode} CLI calculation failed: ${mode_error}")
    endif()
    if(mode STREQUAL "cutoff" AND NOT mode_error MATCHES "Calculating eem")
        message(FATAL_ERROR "cutoff CLI progress output was not emitted: ${mode_error}")
    endif()

    file(READ "${mode_output_prefix}.json" mode_json)
    string(JSON mode_status GET "${mode_json}" results 0 status)
    string(JSON effective_mode GET "${mode_json}" calculation_provenance effective execution mode)
    if(NOT mode_status STREQUAL "success" OR NOT effective_mode STREQUAL "${mode}")
        message(FATAL_ERROR "Unexpected ${mode} CLI result")
    endif()
    file(REMOVE_RECURSE "${mode_output_directory}")
endforeach()

set(warning_output_directory "${CMAKE_CURRENT_BINARY_DIR}/chargefw_cli_warning")
file(REMOVE_RECURSE "${warning_output_directory}")
execute_process(
        COMMAND "${CMAKE_COMMAND}" -E env
                "CHARGEFW_PARAMETER_DIR=${CHARGEFW_PARAMETER_DIR}"
                "${CHARGEFW_CLI}" calculate --method eem --execution full --full-atom-threshold 0
                "${CHARGEFW_INPUT}" "${warning_output_directory}"
        RESULT_VARIABLE warning_result
        ERROR_VARIABLE warning_error
)
if(NOT warning_result EQUAL 0 OR NOT warning_error MATCHES "Warning: ")
    message(FATAL_ERROR "Explicit-full warning was not emitted: ${warning_error}")
endif()
file(REMOVE_RECURSE "${warning_output_directory}")

execute_process(
        COMMAND "${CMAKE_COMMAND}" -E env
                "CHARGEFW_PARAMETER_DIR=${CHARGEFW_PARAMETER_DIR}"
                "${CHARGEFW_CLI}" calculate "${CHARGEFW_INPUT}" "${output_directory}"
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
string(JSON requested_execution GET "${json_output}" calculation_provenance requested execution kind)
if(NOT requested_execution STREQUAL "auto")
    message(FATAL_ERROR "Unexpected requested execution kind: ${requested_execution}")
endif()
string(JSON execution_mode GET "${json_output}" calculation_provenance effective execution mode)
if(NOT execution_mode STREQUAL "full")
    message(FATAL_ERROR "Unexpected effective execution mode: ${execution_mode}")
endif()
string(JSON selected_method GET "${json_output}" calculation_provenance effective method id)
if(selected_method STREQUAL "")
    message(FATAL_ERROR "Selected method provenance is empty")
endif()
string(JSON charge_correction GET "${json_output}" calculation_provenance effective execution charge_correction)
if(NOT charge_correction STREQUAL "none")
    message(FATAL_ERROR "Unexpected full-execution correction policy: ${charge_correction}")
endif()
string(JSON permissive_types GET "${json_output}" calculation_provenance requested classification permissive_types)
if(permissive_types)
    message(FATAL_ERROR "Expected strict type classification by default")
endif()
string(JSON full_atom_threshold GET "${json_output}" calculation_provenance requested resource_policy full_atom_threshold)
if(NOT full_atom_threshold STREQUAL "20000")
    message(FATAL_ERROR "Unexpected default full atom threshold: ${full_atom_threshold}")
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

execute_process(
        COMMAND "${CMAKE_COMMAND}" -E env
                "CHARGEFW_PARAMETER_DIR=${CHARGEFW_PARAMETER_DIR}"
                "${CHARGEFW_CLI}" inspect "${CHARGEFW_INPUT}"
        RESULT_VARIABLE inspect_result
        OUTPUT_VARIABLE inspect_output
)
if(NOT inspect_result EQUAL 0 OR NOT inspect_output MATCHES "records: 1" OR
   NOT inspect_output MATCHES "atoms=3")
    message(FATAL_ERROR "inspect output is incomplete: ${inspect_output}")
endif()

execute_process(
        COMMAND "${CMAKE_COMMAND}" -E env
                "CHARGEFW_PARAMETER_DIR=${CHARGEFW_PARAMETER_DIR}"
                "${CHARGEFW_CLI}" applicability --method formal "${CHARGEFW_INPUT}"
        RESULT_VARIABLE applicability_result
        OUTPUT_VARIABLE applicability_output
)
if(NOT applicability_result EQUAL 0 OR NOT applicability_output MATCHES "applicable method=formal")
    message(FATAL_ERROR "applicability output is incomplete: ${applicability_output}")
endif()

execute_process(
        COMMAND "${CMAKE_COMMAND}" -E env
                "CHARGEFW_PARAMETER_DIR=${CHARGEFW_PARAMETER_DIR}"
                "${CHARGEFW_CLI}" applicability --method smpqeq "${CHARGEFW_INPUT}"
        RESULT_VARIABLE rejected_applicability_result
        OUTPUT_VARIABLE rejected_applicability_output
)
if(NOT rejected_applicability_result EQUAL 0 OR
   NOT rejected_applicability_output MATCHES "rejected method=smpqeq")
    message(FATAL_ERROR "rejected applicability output has incorrect method identity: ${rejected_applicability_output}")
endif()

execute_process(COMMAND "${CHARGEFW_CLI}" methods OUTPUT_VARIABLE methods_output)
if(NOT methods_output MATCHES "formal" OR NOT methods_output MATCHES "minimum>=1" OR
   NOT methods_output MATCHES "minimum>0")
    message(FATAL_ERROR "methods output does not list formal")
endif()

execute_process(
        COMMAND "${CMAKE_COMMAND}" -E env
                "CHARGEFW_PARAMETER_DIR=${CHARGEFW_PARAMETER_DIR}"
                "${CHARGEFW_CLI}" parameters qeq
        OUTPUT_VARIABLE parameters_output
)
if(NOT parameters_output MATCHES "QEq_original" OR parameters_output MATCHES "EEM_original")
    message(FATAL_ERROR "parameter filtering is incorrect: ${parameters_output}")
endif()
