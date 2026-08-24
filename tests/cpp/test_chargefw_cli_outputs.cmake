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
    if(mode STREQUAL "cutoff")
        string(ASCII 27 erase_escape)
        set(clear_line "\r${erase_escape}[2K")
        string(REPLACE "${clear_line}" "" progress_without_clear_sequences "${mode_error}")
        if(progress_without_clear_sequences MATCHES "\r")
            message(FATAL_ERROR "progress output did not erase a previous rendered line: ${mode_error}")
        endif()
        string(FIND "${mode_error}" "${clear_line}\n" final_clear_position)
        if(NOT mode_error MATCHES "Fragments" OR NOT mode_error MATCHES "Targets" OR
           final_clear_position EQUAL -1)
            message(FATAL_ERROR "progress output did not clear its final rendered line: ${mode_error}")
        endif()
    endif()

    file(READ "${mode_output_prefix}.json" mode_json)
    string(JSON mode_status GET "${mode_json}" results 0 status)
    string(JSON effective_mode GET "${mode_json}" calculation_provenance effective execution mode)
    if(NOT mode_status STREQUAL "success" OR NOT effective_mode STREQUAL "${mode}")
        message(FATAL_ERROR "Unexpected ${mode} CLI result")
    endif()
    file(REMOVE_RECURSE "${mode_output_directory}")
endforeach()

set(no_plan_output_directory "${CMAKE_CURRENT_BINARY_DIR}/chargefw_cli_no_plan")
set(no_plan_output_prefix "${no_plan_output_directory}/water.chargefw")
file(REMOVE_RECURSE "${no_plan_output_directory}")
execute_process(
        COMMAND "${CMAKE_COMMAND}" -E env
                "CHARGEFW_PARAMETER_DIR=${CHARGEFW_PARAMETER_DIR}"
                "${CHARGEFW_CLI}" calculate --method smpqeq "${CHARGEFW_INPUT}"
                "${no_plan_output_directory}"
        RESULT_VARIABLE no_plan_result
        OUTPUT_VARIABLE no_plan_output
        ERROR_VARIABLE no_plan_error
)
if(NOT no_plan_result EQUAL 3)
    message(FATAL_ERROR "no-plan CLI exit status was ${no_plan_result}: ${no_plan_error}")
endif()
if(NOT EXISTS "${no_plan_output_prefix}.json")
    message(FATAL_ERROR "no-plan CLI JSON output was not written to the output directory")
endif()
foreach(extension IN ITEMS sdf mol2 cif)
    if(EXISTS "${no_plan_output_prefix}.${extension}")
        message(FATAL_ERROR "no-plan CLI wrote unexpected molecular output: ${extension}")
    endif()
endforeach()
file(READ "${no_plan_output_prefix}.json" no_plan_json)
string(JSON no_plan_status GET "${no_plan_json}" status)
string(JSON no_plan_record_status GET "${no_plan_json}" results 0 status)
string(JSON no_plan_diagnostic GET "${no_plan_json}" diagnostics 0 code)
if(NOT no_plan_status STREQUAL "no_executable_plan" OR
   NOT no_plan_record_status STREQUAL "no_executable_plan" OR
   NOT no_plan_diagnostic STREQUAL "no_executable_plan")
    message(FATAL_ERROR "no-plan CLI JSON does not report a structured no-plan result")
endif()
file(REMOVE_RECURSE "${no_plan_output_directory}")

set(invalid_output_directory "${CMAKE_CURRENT_BINARY_DIR}/chargefw_cli_invalid")
set(invalid_output_prefix "${invalid_output_directory}/water.chargefw")
file(REMOVE_RECURSE "${invalid_output_directory}")
execute_process(
        COMMAND "${CMAKE_COMMAND}" -E env
                "CHARGEFW_PARAMETER_DIR=${CHARGEFW_PARAMETER_DIR}"
                "${CHARGEFW_CLI}" calculate --threads 18446744073709551615 "${CHARGEFW_INPUT}"
                "${invalid_output_directory}"
        RESULT_VARIABLE invalid_result
        ERROR_VARIABLE invalid_error
)
if(NOT invalid_result EQUAL 2)
    message(FATAL_ERROR "invalid-request CLI exit status was ${invalid_result}: ${invalid_error}")
endif()
if(NOT EXISTS "${invalid_output_prefix}.json")
    message(FATAL_ERROR "invalid-request CLI JSON output was not written to the output directory")
endif()
file(READ "${invalid_output_prefix}.json" invalid_json)
string(JSON invalid_status GET "${invalid_json}" status)
string(JSON invalid_diagnostic GET "${invalid_json}" diagnostics 0 code)
if(NOT invalid_status STREQUAL "invalid_input_or_request" OR
   NOT invalid_diagnostic STREQUAL "invalid_input_or_request")
    message(FATAL_ERROR "invalid-request CLI JSON does not report a structured invalid result")
endif()
file(REMOVE_RECURSE "${invalid_output_directory}")

set(warning_output_directory "${CMAKE_CURRENT_BINARY_DIR}/chargefw_cli_warning")
file(REMOVE_RECURSE "${warning_output_directory}")
execute_process(
        COMMAND "${CMAKE_COMMAND}" -E env
                "CHARGEFW_PARAMETER_DIR=${CHARGEFW_PARAMETER_DIR}"
                "${CHARGEFW_CLI}" calculate --method eem --execution full --cutoff-atom-threshold 0
                --cover-atom-threshold 0
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
string(JSON cutoff_atom_threshold GET "${json_output}" calculation_provenance requested resource_policy cutoff_atom_threshold)
if(NOT cutoff_atom_threshold STREQUAL "20000")
    message(FATAL_ERROR "Unexpected default cutoff atom threshold: ${cutoff_atom_threshold}")
endif()
string(JSON cover_atom_threshold GET "${json_output}" calculation_provenance requested resource_policy cover_atom_threshold)
if(NOT cover_atom_threshold STREQUAL "80000")
    message(FATAL_ERROR "Unexpected default cover atom threshold: ${cover_atom_threshold}")
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
